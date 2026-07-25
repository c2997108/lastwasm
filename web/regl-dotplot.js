const FORWARD_COLOR = [22 / 255, 100 / 255, 217 / 255, 1];
const REVERSE_COLOR = [201 / 255, 42 / 255, 42 / 255, 1];
const PLOT_MARGIN = { left: 68, right: 18, top: 22, bottom: 54 };
const PLOTLY_VISIBLE_ALIGNMENT_TARGET = 100000;

const VISUAL_LINE_VERTEX_SHADER = `
precision highp float;
attribute vec2 position;
attribute float strand;
attribute float alignmentId;
uniform vec2 viewCenter;
uniform vec2 plotSize;
uniform float pixelsPerUnit;
uniform float pointSize;
varying float vStrand;
varying float vAlignmentId;

void main() {
  vec2 clipPosition = (position - viewCenter) * pixelsPerUnit * 2.0 / plotSize;
  gl_Position = vec4(clipPosition, 0.0, 1.0);
  gl_PointSize = pointSize;
  vStrand = strand;
  vAlignmentId = alignmentId;
}
`;

const PICK_VERTEX_SHADER = `
precision highp float;
attribute vec2 position;
attribute float alignmentId;
uniform vec2 viewCenter;
uniform vec2 plotSize;
uniform float pixelsPerUnit;
uniform float pointSize;
varying float vAlignmentId;
void main() {
  vec2 clipPosition = (position - viewCenter) * pixelsPerUnit * 2.0 / plotSize;
  gl_Position = vec4(clipPosition, 0.0, 1.0);
  gl_PointSize = pointSize;
  vAlignmentId = alignmentId;
}
`;

const VISUAL_FRAGMENT_SHADER = `
precision mediump float;
uniform vec4 forwardColor;
uniform vec4 reverseColor;
uniform float representativeEvery;
uniform float backgroundOpacity;
varying float vStrand;
varying float vAlignmentId;
void main() {
  vec4 color = mix(forwardColor, reverseColor, step(0.5, vStrand));
  float representative = 1.0 - step(0.5, mod(floor(vAlignmentId + 0.5), representativeEvery));
  float opacity = mix(backgroundOpacity, 1.0, representative);
  gl_FragColor = vec4(color.rgb, opacity);
}
`;

const PICK_FRAGMENT_SHADER = `
precision highp float;
varying float vAlignmentId;
void main() {
  float value = floor(vAlignmentId + 1.5);
  float red = mod(value, 256.0);
  float green = mod(floor(value / 256.0), 256.0);
  float blue = mod(floor(value / 65536.0), 256.0);
  gl_FragColor = vec4(red, green, blue, 255.0) / 255.0;
}
`;

export async function renderReglDotPlot(container, plotData) {
  const renderer = new ReglDotPlot(container, plotData);
  const timing = await renderer.initialize();
  return { renderer, timing };
}

class ReglDotPlot {
  constructor(container, plotData) {
    this.container = container;
    this.data = plotData;
    this.regl = null;
    this.gl = null;
    this.resources = [];
    this.pickingFramebuffer = null;
    this.pickingDirty = true;
    this.destroyed = false;
    this.drawRequest = 0;
    this.pickRequest = 0;
    this.pendingPointer = null;
    this.dragState = null;
    this.viewModified = false;
    this.dpr = 1;
    this.plotRect = { left: 0, top: 0, width: 1, height: 1 };
    this.view = { centerX: plotData.maxX / 2, centerY: plotData.maxY / 2, scale: 1 };
    this.representativeEvery = Math.max(1, Math.ceil(plotData.renderedCount / PLOTLY_VISIBLE_ALIGNMENT_TARGET));
    this.backgroundOpacity = this.representativeEvery === 1
      ? 1
      : Math.min(0.08, 0.5 / this.representativeEvery);
  }

  async initialize() {
    const totalStartedAt = performance.now();
    this.buildDom();
    this.createContext();
    this.createCommands();
    const uploadStartedAt = performance.now();
    this.uploadChunks();
    const uploadMs = performance.now() - uploadStartedAt;
    this.resize(true);
    const drawStartedAt = performance.now();
    this.drawNow();
    const drawSubmitMs = performance.now() - drawStartedAt;
    const gpuStartedAt = performance.now();
    await waitForGpu(this.gl);
    const gpuCompleteMs = performance.now() - gpuStartedAt;
    const paintStartedAt = performance.now();
    await afterTwoFrames();
    const firstPaintMs = performance.now() - paintStartedAt;
    this.installInteraction();
    this.resizeObserver = new ResizeObserver(() => this.resize(false));
    this.resizeObserver.observe(this.container);
    return {
      uploadMs,
      drawSubmitMs,
      gpuCompleteMs,
      firstPaintMs,
      totalMs: performance.now() - totalStartedAt,
    };
  }

  buildDom() {
    this.container.replaceChildren();
    this.stage = document.createElement('div');
    this.stage.className = 'dotplot-stage';
    this.stage.setAttribute('aria-label', 'Alignment dot plot');
    this.glCanvas = document.createElement('canvas');
    this.glCanvas.className = 'dotplot-gl';
    this.axisCanvas = document.createElement('canvas');
    this.axisCanvas.className = 'dotplot-axis';
    this.tooltip = document.createElement('div');
    this.tooltip.className = 'dotplot-tooltip';
    this.tooltip.hidden = true;
    this.fitButton = document.createElement('button');
    this.fitButton.className = 'dotplot-fit';
    this.fitButton.type = 'button';
    this.fitButton.textContent = 'Fit';
    this.fitButton.title = 'Fit all alignments';
    this.fitButton.addEventListener('pointerdown', event => event.stopPropagation());
    this.fitButton.addEventListener('click', () => this.fitView());
    this.stage.append(this.glCanvas, this.axisCanvas, this.tooltip, this.fitButton);
    this.container.append(this.stage);
  }

  createContext() {
    if (typeof window.createREGL !== 'function') throw new Error('regl failed to load.');
    this.gl = this.glCanvas.getContext('webgl2', {
      alpha: false,
      antialias: true,
      depth: false,
      preserveDrawingBuffer: false,
    });
    if (!this.gl) throw new Error('WebGL2 is required for the full-resolution dot plot.');
    this.gl.disable(this.gl.DITHER);
    this.regl = window.createREGL({ gl: this.gl });
  }

  createCommands() {
    const vertexAttributes = {
      position: {
        buffer: this.regl.prop('positionBuffer'),
        size: 2,
      },
      strand: {
        buffer: this.regl.prop('strandBuffer'),
        size: 1,
      },
      alignmentId: {
        buffer: this.regl.prop('idBuffer'),
        size: 1,
      },
    };
    const uniforms = {
      viewCenter: this.regl.prop('viewCenter'),
      plotSize: this.regl.prop('plotSize'),
      pixelsPerUnit: this.regl.prop('pixelsPerUnit'),
      pointSize: this.regl.prop('pointSize'),
    };
    const viewport = {
      x: this.regl.prop('viewportX'),
      y: this.regl.prop('viewportY'),
      width: this.regl.prop('viewportWidth'),
      height: this.regl.prop('viewportHeight'),
    };
    const visualAttributes = {
      ...vertexAttributes,
    };
    const visualUniforms = {
      viewCenter: uniforms.viewCenter,
      plotSize: uniforms.plotSize,
      pixelsPerUnit: uniforms.pixelsPerUnit,
      pointSize: this.regl.prop('pointSize'),
      forwardColor: FORWARD_COLOR,
      reverseColor: REVERSE_COLOR,
      representativeEvery: this.regl.prop('representativeEvery'),
      backgroundOpacity: this.regl.prop('backgroundOpacity'),
    };
    const visualState = {
      attributes: visualAttributes,
      uniforms: visualUniforms,
      viewport,
      count: this.regl.prop('vertexCount'),
      depth: { enable: false },
      blend: {
        enable: true,
        func: { srcRGB: 'src alpha', dstRGB: 'one minus src alpha', srcAlpha: 'one', dstAlpha: 'one minus src alpha' },
      },
      cull: { enable: false },
    };
    this.drawVisualLines = this.regl({
      vert: VISUAL_LINE_VERTEX_SHADER,
      frag: VISUAL_FRAGMENT_SHADER,
      ...visualState,
      primitive: 'lines',
      lineWidth: 1,
    });
    this.drawVisualPoints = this.regl({
      vert: VISUAL_LINE_VERTEX_SHADER,
      frag: VISUAL_FRAGMENT_SHADER,
      ...visualState,
      primitive: 'points',
    });
    const pickingState = {
      vert: PICK_VERTEX_SHADER,
      frag: PICK_FRAGMENT_SHADER,
      attributes: {
        position: vertexAttributes.position,
        alignmentId: vertexAttributes.alignmentId,
      },
      uniforms: {
        viewCenter: uniforms.viewCenter,
        plotSize: uniforms.plotSize,
        pixelsPerUnit: uniforms.pixelsPerUnit,
        pointSize: uniforms.pointSize,
      },
      viewport,
      framebuffer: this.regl.prop('framebuffer'),
      count: this.regl.prop('vertexCount'),
      depth: { enable: false },
      blend: { enable: false },
      cull: { enable: false },
    };
    this.drawPickingLines = this.regl({
      ...pickingState,
      primitive: 'lines',
      lineWidth: 1,
    });
    this.drawPickingPoints = this.regl({
      ...pickingState,
      primitive: 'points',
    });
  }

  uploadChunks() {
    for (const chunk of this.data.chunks) {
      const arrays = chunk.arrays;
      const refNameIds = new (chunk.nameIdBytes === 2 ? Uint16Array : Uint32Array)(arrays.refNameIds);
      const positionBuffer = this.regl.buffer({ data: new Float32Array(arrays.positions), usage: 'static' });
      const strandBuffer = this.regl.buffer({ data: new Uint8Array(arrays.strands), type: 'uint8', usage: 'static' });
      const idBuffer = this.regl.buffer({ data: new Float32Array(arrays.alignmentIds), usage: 'static' });
      const NameIdArray = chunk.nameIdBytes === 2 ? Uint16Array : Uint32Array;
      this.resources.push({
        index: chunk.index,
        alignmentStart: chunk.alignmentStart,
        alignmentCount: chunk.alignmentCount,
        positionBuffer,
        strandBuffer,
        idBuffer,
        metadata: {
          scores: new Float64Array(arrays.scores),
          refNameIds,
          qryNameIds: new NameIdArray(arrays.qryNameIds),
          refStarts: new Uint32Array(arrays.refStarts),
          refLengths: new Uint32Array(arrays.refLengths),
          qryStarts: new Uint32Array(arrays.qryStarts),
          qryLengths: new Uint32Array(arrays.qryLengths),
          strands: new Uint8Array(arrays.strands),
        },
      });
      arrays.positions = null;
      arrays.alignmentIds = null;
    }
    this.resources.sort((left, right) => left.alignmentStart - right.alignmentStart);
  }

  resize(refit) {
    if (this.destroyed) return;
    const bounds = this.container.getBoundingClientRect();
    const width = Math.max(320, Math.round(bounds.width));
    const height = Math.max(320, Math.round(bounds.height));
    this.dpr = Math.max(1, Math.min(1.5, window.devicePixelRatio || 1));
    const pixelWidth = Math.round(width * this.dpr);
    const pixelHeight = Math.round(height * this.dpr);
    if (this.glCanvas.width !== pixelWidth || this.glCanvas.height !== pixelHeight) {
      this.glCanvas.width = pixelWidth;
      this.glCanvas.height = pixelHeight;
      this.axisCanvas.width = pixelWidth;
      this.axisCanvas.height = pixelHeight;
      this.glCanvas.style.width = `${width}px`;
      this.glCanvas.style.height = `${height}px`;
      this.axisCanvas.style.width = `${width}px`;
      this.axisCanvas.style.height = `${height}px`;
      this.destroyPickingFramebuffer();
    }
    const margin = width < 600
      ? { left: 54, right: 12, top: 18, bottom: 48 }
      : PLOT_MARGIN;
    this.plotRect = {
      left: margin.left,
      top: margin.top,
      width: Math.max(1, width - margin.left - margin.right),
      height: Math.max(1, height - margin.top - margin.bottom),
    };
    if (refit || !this.viewModified) this.fitView(false);
    this.requestDraw();
  }

  fitView(draw = true) {
    const maxX = Math.max(1, this.data.maxX);
    const maxY = Math.max(1, this.data.maxY);
    this.view.centerX = maxX / 2;
    this.view.centerY = maxY / 2;
    this.view.scale = Math.max(1e-12, Math.min(this.plotRect.width / maxX, this.plotRect.height / maxY) * 0.96);
    this.fitScale = this.view.scale;
    this.viewModified = false;
    this.stage.dataset.viewMode = 'fit';
    this.hideTooltip();
    if (draw) this.requestDraw();
  }

  drawProps(resource, framebuffer = null) {
    const rect = this.plotRect;
    const zoomOpacity = this.backgroundOpacity * Math.max(1, this.view.scale / Math.max(this.fitScale || this.view.scale, 1e-12));
    return {
      positionBuffer: resource.positionBuffer,
      strandBuffer: resource.strandBuffer,
      idBuffer: resource.idBuffer,
      vertexCount: resource.alignmentCount * 2,
      viewCenter: [this.view.centerX, this.view.centerY],
      plotSize: [rect.width, rect.height],
      pixelsPerUnit: this.view.scale,
      pointSize: Math.max(2, 1.5 * this.dpr),
      representativeEvery: this.representativeEvery,
      backgroundOpacity: Math.min(1, zoomOpacity),
      viewportX: Math.round(rect.left * this.dpr),
      viewportY: Math.round((this.glCanvas.height / this.dpr - rect.top - rect.height) * this.dpr),
      viewportWidth: Math.max(1, Math.round(rect.width * this.dpr)),
      viewportHeight: Math.max(1, Math.round(rect.height * this.dpr)),
      framebuffer,
    };
  }

  drawNow() {
    if (this.destroyed) return;
    this.regl.clear({ color: [1, 1, 1, 1], depth: 1 });
    for (const resource of this.resources) {
      const props = this.drawProps(resource);
      this.drawVisualLines(props);
      this.drawVisualPoints(props);
    }
    this.drawAxes();
    this.destroyPickingFramebuffer();
  }

  requestDraw() {
    if (this.drawRequest || this.destroyed) return;
    this.drawRequest = requestAnimationFrame(() => {
      this.drawRequest = 0;
      this.drawNow();
    });
  }

  drawAxes() {
    const canvas = this.axisCanvas;
    const ctx = canvas.getContext('2d');
    const width = canvas.width / this.dpr;
    const height = canvas.height / this.dpr;
    ctx.setTransform(this.dpr, 0, 0, this.dpr, 0, 0);
    ctx.clearRect(0, 0, width, height);
    const rect = this.plotRect;
    const range = this.visibleRange();
    const xTicks = makeTicks(range.minX, range.maxX, 6);
    const yTicks = makeTicks(range.minY, range.maxY, 6);
    const xTickLabels = formatAxisNumbers(xTicks);
    const yTickLabels = formatAxisNumbers(yTicks);
    ctx.save();
    ctx.beginPath();
    ctx.rect(rect.left, rect.top, rect.width, rect.height);
    ctx.clip();
    this.drawSequenceBoundaries(ctx, true, range.minX, range.maxX);
    this.drawSequenceBoundaries(ctx, false, range.minY, range.maxY);
    ctx.restore();

    ctx.strokeStyle = '#8b95a5';
    ctx.lineWidth = 1;
    ctx.setLineDash([]);
    ctx.strokeRect(rect.left + 0.5, rect.top + 0.5, rect.width - 1, rect.height - 1);
    ctx.beginPath();
    for (const tick of xTicks) {
      const x = this.dataToScreenX(tick);
      if (x < rect.left || x > rect.left + rect.width) continue;
      ctx.moveTo(x, rect.top + rect.height);
      ctx.lineTo(x, rect.top + rect.height + 4);
    }
    for (const tick of yTicks) {
      const y = this.dataToScreenY(tick);
      if (y < rect.top || y > rect.top + rect.height) continue;
      ctx.moveTo(rect.left - 4, y);
      ctx.lineTo(rect.left, y);
    }
    ctx.stroke();
    ctx.fillStyle = '#4d596a';
    ctx.font = '11px system-ui, sans-serif';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'top';
    for (let index = 0; index < xTicks.length; index++) {
      const tick = xTicks[index];
      const x = this.dataToScreenX(tick);
      if (x < rect.left || x > rect.left + rect.width) continue;
      ctx.fillText(xTickLabels[index], x, rect.top + rect.height + 6);
    }
    ctx.textAlign = 'right';
    ctx.textBaseline = 'middle';
    for (let index = 0; index < yTicks.length; index++) {
      const tick = yTicks[index];
      const y = this.dataToScreenY(tick);
      if (y < rect.top || y > rect.top + rect.height) continue;
      ctx.fillText(yTickLabels[index], rect.left - 7, y);
    }
    ctx.fillStyle = '#222832';
    ctx.font = '12px system-ui, sans-serif';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'bottom';
    ctx.fillText('Reference', rect.left + rect.width / 2, height - 2);
    ctx.save();
    ctx.translate(13, rect.top + rect.height / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.fillText('Query', 0, 0);
    ctx.restore();
    this.drawLegend(ctx, rect);
  }

  drawSequenceBoundaries(ctx, reference, min, max) {
    ctx.save();
    const orders = reference ? this.data.refOrders : this.data.qryOrders;
    const offsets = reference ? this.data.refOffsets : this.data.qryOffsets;
    const lengths = reference ? this.data.refLenByName : this.data.qryLenByName;
    ctx.strokeStyle = 'rgba(0, 0, 0, 0.14)';
    ctx.setLineDash([5, 4]);
    ctx.fillStyle = '#596579';
    ctx.font = '10px system-ui, sans-serif';
    for (let index = 0; index < orders.length; index++) {
      const name = orders[index];
      const start = offsets.get(name) || 0;
      const end = start + (lengths.get(name) || 0);
      if (end < min || start > max) continue;
      if (start > 0) {
        const position = reference ? this.dataToScreenX(start) : this.dataToScreenY(start);
        ctx.beginPath();
        if (reference) {
          ctx.moveTo(position, this.plotRect.top);
          ctx.lineTo(position, this.plotRect.top + this.plotRect.height);
        } else {
          ctx.moveTo(this.plotRect.left, position);
          ctx.lineTo(this.plotRect.left + this.plotRect.width, position);
        }
        ctx.stroke();
      }
      const spanPixels = (end - start) * this.view.scale;
      if (spanPixels < 58) continue;
      const midpoint = (start + end) / 2;
      if (reference) {
        ctx.textAlign = 'center';
        ctx.textBaseline = 'top';
        ctx.fillText(name, this.dataToScreenX(midpoint), this.plotRect.top + 4, Math.max(20, spanPixels - 6));
      } else {
        const y = this.dataToScreenY(midpoint);
        ctx.save();
        ctx.translate(this.plotRect.left + 4, y);
        ctx.rotate(-Math.PI / 2);
        ctx.textAlign = 'center';
        ctx.textBaseline = 'top';
        ctx.fillText(name, 0, 0, Math.max(20, spanPixels - 6));
        ctx.restore();
      }
    }
    ctx.restore();
  }

  drawLegend(ctx, rect) {
    const x = rect.left + rect.width - 112;
    const y = rect.top + 9;
    ctx.fillStyle = 'rgba(255,255,255,0.88)';
    ctx.fillRect(x - 7, y - 5, 112, 38);
    ctx.font = '11px system-ui, sans-serif';
    ctx.textAlign = 'left';
    ctx.textBaseline = 'middle';
    ctx.strokeStyle = '#1664d9';
    ctx.lineWidth = 2;
    ctx.beginPath(); ctx.moveTo(x, y + 3); ctx.lineTo(x + 18, y + 3); ctx.stroke();
    ctx.fillStyle = '#333b47'; ctx.fillText('Forward', x + 24, y + 3);
    ctx.strokeStyle = '#c92a2a';
    ctx.beginPath(); ctx.moveTo(x, y + 21); ctx.lineTo(x + 18, y + 21); ctx.stroke();
    ctx.fillText('Reverse', x + 24, y + 21);
  }

  visibleRange() {
    const halfX = this.plotRect.width / (2 * this.view.scale);
    const halfY = this.plotRect.height / (2 * this.view.scale);
    return {
      minX: this.view.centerX - halfX,
      maxX: this.view.centerX + halfX,
      minY: this.view.centerY - halfY,
      maxY: this.view.centerY + halfY,
    };
  }

  dataToScreenX(value) {
    return this.plotRect.left + this.plotRect.width / 2 + (value - this.view.centerX) * this.view.scale;
  }

  dataToScreenY(value) {
    return this.plotRect.top + this.plotRect.height / 2 - (value - this.view.centerY) * this.view.scale;
  }

  screenToData(x, y) {
    return {
      x: this.view.centerX + (x - this.plotRect.left - this.plotRect.width / 2) / this.view.scale,
      y: this.view.centerY - (y - this.plotRect.top - this.plotRect.height / 2) / this.view.scale,
    };
  }

  installInteraction() {
    this.onPointerDown = event => {
      const point = this.pointerPosition(event);
      if (!this.containsPlotPoint(point)) return;
      this.dragState = {
        pointerId: event.pointerId,
        startX: point.x,
        startY: point.y,
        centerX: this.view.centerX,
        centerY: this.view.centerY,
      };
      this.stage.setPointerCapture(event.pointerId);
      this.stage.classList.add('is-panning');
      this.hideTooltip();
    };
    this.onPointerMove = event => {
      const point = this.pointerPosition(event);
      if (this.dragState?.pointerId === event.pointerId) {
        this.view.centerX = this.dragState.centerX - (point.x - this.dragState.startX) / this.view.scale;
        this.view.centerY = this.dragState.centerY + (point.y - this.dragState.startY) / this.view.scale;
        this.viewModified = true;
        this.stage.dataset.viewMode = 'custom';
        this.requestDraw();
        return;
      }
      if (!this.containsPlotPoint(point)) {
        this.hideTooltip();
        return;
      }
      this.pendingPointer = point;
      if (!this.pickRequest) {
        this.pickRequest = requestAnimationFrame(() => {
          this.pickRequest = 0;
          this.pickAt(this.pendingPointer);
        });
      }
    };
    this.onPointerUp = event => {
      if (this.dragState?.pointerId !== event.pointerId) return;
      this.dragState = null;
      this.stage.classList.remove('is-panning');
      this.stage.releasePointerCapture(event.pointerId);
    };
    this.onWheel = event => {
      const point = this.pointerPosition(event);
      if (!this.containsPlotPoint(point)) return;
      event.preventDefault();
      const before = this.screenToData(point.x, point.y);
      const factor = Math.exp(-event.deltaY * 0.0015);
      const minScale = Math.max(1e-12, this.fitScale * 0.1);
      const maxScale = this.fitScale * 1e7;
      this.view.scale = Math.max(minScale, Math.min(maxScale, this.view.scale * factor));
      this.view.centerX = before.x - (point.x - this.plotRect.left - this.plotRect.width / 2) / this.view.scale;
      this.view.centerY = before.y + (point.y - this.plotRect.top - this.plotRect.height / 2) / this.view.scale;
      this.viewModified = true;
      this.stage.dataset.viewMode = 'custom';
      this.hideTooltip();
      this.requestDraw();
    };
    this.onDoubleClick = event => {
      if (this.containsPlotPoint(this.pointerPosition(event))) this.fitView();
    };
    this.stage.addEventListener('pointerdown', this.onPointerDown);
    this.stage.addEventListener('pointermove', this.onPointerMove);
    this.stage.addEventListener('pointerup', this.onPointerUp);
    this.stage.addEventListener('pointercancel', this.onPointerUp);
    this.stage.addEventListener('wheel', this.onWheel, { passive: false });
    this.stage.addEventListener('dblclick', this.onDoubleClick);
  }

  pointerPosition(event) {
    const bounds = this.stage.getBoundingClientRect();
    return { x: event.clientX - bounds.left, y: event.clientY - bounds.top };
  }

  containsPlotPoint(point) {
    const rect = this.plotRect;
    return point && point.x >= rect.left && point.x <= rect.left + rect.width
      && point.y >= rect.top && point.y <= rect.top + rect.height;
  }

  ensurePickingBuffer() {
    if (!this.pickingFramebuffer) {
      this.pickingFramebuffer = this.regl.framebuffer({
        color: this.regl.texture({
          width: this.glCanvas.width,
          height: this.glCanvas.height,
          format: 'rgba',
          type: 'uint8',
          min: 'nearest',
          mag: 'nearest',
        }),
        depth: false,
        stencil: false,
      });
      this.pickingDirty = true;
    }
    if (!this.pickingDirty) return;
    this.regl.clear({ framebuffer: this.pickingFramebuffer, color: [0, 0, 0, 0] });
    for (const resource of this.resources) {
      const props = this.drawProps(resource, this.pickingFramebuffer);
      this.drawPickingLines(props);
      this.drawPickingPoints(props);
    }
    this.pickingDirty = false;
  }

  pickAt(point) {
    if (this.destroyed || !point || !this.containsPlotPoint(point)) return;
    this.ensurePickingBuffer();
    const radius = Math.max(2, Math.round(7 * this.dpr));
    const centerX = Math.floor(point.x * this.dpr);
    const centerY = this.glCanvas.height - 1 - Math.floor(point.y * this.dpr);
    const x = Math.max(0, centerX - radius);
    const y = Math.max(0, centerY - radius);
    const width = Math.min(this.glCanvas.width - x, radius * 2 + 1);
    const height = Math.min(this.glCanvas.height - y, radius * 2 + 1);
    const pixels = new Uint8Array(width * height * 4);
    this.regl.read({ framebuffer: this.pickingFramebuffer, x, y, width, height, data: pixels });
    const encoded = nearestEncodedPixel(pixels, width, height, centerX - x, centerY - y);
    if (encoded === 0) {
      this.hideTooltip();
      return;
    }
    const alignmentId = encoded - 1;
    const detail = this.alignmentDetail(alignmentId);
    if (!detail) {
      this.hideTooltip();
      return;
    }
    this.tooltip.textContent = [
      `Score: ${detail.score}`,
      `Ref: ${detail.refName}:${detail.refStart}-${detail.refEnd}`,
      `Qry: ${detail.qryName}:${detail.qryStart}-${detail.qryEnd}`,
      `Length: ${detail.length}`,
    ].join('\n');
    this.tooltip.hidden = false;
    const stageWidth = this.stage.clientWidth;
    const stageHeight = this.stage.clientHeight;
    const tooltipWidth = this.tooltip.offsetWidth;
    const tooltipHeight = this.tooltip.offsetHeight;
    this.tooltip.style.left = `${Math.max(6, Math.min(stageWidth - tooltipWidth - 6, point.x + 14))}px`;
    this.tooltip.style.top = `${Math.max(6, Math.min(stageHeight - tooltipHeight - 6, point.y + 14))}px`;
  }

  alignmentDetail(alignmentId) {
    let low = 0;
    let high = this.resources.length - 1;
    while (low <= high) {
      const middle = (low + high) >> 1;
      const resource = this.resources[middle];
      if (alignmentId < resource.alignmentStart) {
        high = middle - 1;
      } else if (alignmentId >= resource.alignmentStart + resource.alignmentCount) {
        low = middle + 1;
      } else {
        const index = alignmentId - resource.alignmentStart;
        const metadata = resource.metadata;
        const refId = metadata.refNameIds[index];
        const qryId = metadata.qryNameIds[index];
        const refName = this.data.refOrders[refId] || '';
        const qryName = this.data.qryOrders[qryId] || '';
        const refStart = metadata.refStarts[index] + 1;
        const refLength = metadata.refLengths[index];
        const qryRawStart = metadata.qryStarts[index];
        const qryLength = metadata.qryLengths[index];
        const reverse = metadata.strands[index * 2] === 1;
        const qryTotal = this.data.qryLenByName.get(qryName) || 0;
        return {
          score: metadata.scores[index],
          refName,
          refStart,
          refEnd: refStart + refLength - 1,
          qryName,
          qryStart: reverse ? qryTotal - qryRawStart : qryRawStart + 1,
          qryEnd: reverse ? qryTotal - qryRawStart - qryLength + 1 : qryRawStart + qryLength,
          length: Math.min(refLength, qryLength),
        };
      }
    }
    return null;
  }

  hideTooltip() {
    if (this.tooltip) this.tooltip.hidden = true;
  }

  destroyPickingFramebuffer() {
    this.pickingFramebuffer?.destroy();
    this.pickingFramebuffer = null;
    this.pickingDirty = true;
  }

  destroy() {
    if (this.destroyed) return;
    this.destroyed = true;
    this.resizeObserver?.disconnect();
    if (this.drawRequest) cancelAnimationFrame(this.drawRequest);
    if (this.pickRequest) cancelAnimationFrame(this.pickRequest);
    this.stage?.removeEventListener('pointerdown', this.onPointerDown);
    this.stage?.removeEventListener('pointermove', this.onPointerMove);
    this.stage?.removeEventListener('pointerup', this.onPointerUp);
    this.stage?.removeEventListener('pointercancel', this.onPointerUp);
    this.stage?.removeEventListener('wheel', this.onWheel);
    this.stage?.removeEventListener('dblclick', this.onDoubleClick);
    this.destroyPickingFramebuffer();
    for (const resource of this.resources) {
      resource.positionBuffer.destroy();
      resource.strandBuffer.destroy();
      resource.idBuffer.destroy();
    }
    this.regl?.destroy();
    this.container.replaceChildren();
  }
}

function nearestEncodedPixel(pixels, width, height, centerX, centerY) {
  let nearest = 0;
  let nearestDistance = Infinity;
  for (let row = 0; row < height; row++) {
    for (let column = 0; column < width; column++) {
      const index = (row * width + column) * 4;
      if (pixels[index + 3] === 0) continue;
      const distance = (column - centerX) ** 2 + (row - centerY) ** 2;
      if (distance >= nearestDistance) continue;
      nearest = pixels[index] + pixels[index + 1] * 256 + pixels[index + 2] * 65536;
      nearestDistance = distance;
    }
  }
  return nearest;
}

function makeTicks(min, max, count) {
  if (!Number.isFinite(min) || !Number.isFinite(max) || max <= min) return [];
  const rough = (max - min) / Math.max(1, count);
  const magnitude = 10 ** Math.floor(Math.log10(rough));
  const normalized = rough / magnitude;
  const step = (normalized >= 5 ? 5 : normalized >= 2 ? 2 : 1) * magnitude;
  const first = Math.ceil(min / step) * step;
  const ticks = [];
  for (let value = first; value <= max && ticks.length < 20; value += step) ticks.push(value);
  return ticks;
}

export function formatAxisNumbers(values) {
  if (!values.length) return [];
  const maximum = values.reduce((result, value) => Math.max(result, Math.abs(value)), 0);
  const [divisor, suffix] = maximum >= 1e9
    ? [1e9, 'G']
    : maximum >= 1e6
      ? [1e6, 'M']
      : maximum >= 1e3
        ? [1e3, 'k']
        : [1, ''];
  let step = Infinity;
  for (let index = 1; index < values.length; index++) {
    const difference = Math.abs(values[index] - values[index - 1]);
    if (difference > 0) step = Math.min(step, difference);
  }
  const scaledStep = Number.isFinite(step) ? step / divisor : 1;
  const requiredDecimals = scaledStep > 0
    ? Math.max(0, Math.ceil(-Math.log10(scaledStep) - 1e-10))
    : 0;
  const decimals = Math.min(9, suffix ? Math.max(1, requiredDecimals) : requiredDecimals);
  const zeroTolerance = Number.isFinite(step) ? step * 1e-9 : 0;
  return values.map(value => {
    const normalized = Math.abs(value) < zeroTolerance ? 0 : value / divisor;
    return `${normalized.toFixed(decimals)}${suffix}`;
  });
}

function waitForGpu(gl) {
  if (typeof gl.fenceSync !== 'function') return Promise.resolve();
  const sync = gl.fenceSync(gl.SYNC_GPU_COMMANDS_COMPLETE, 0);
  gl.flush();
  return new Promise(resolve => {
    const poll = () => {
      const status = gl.clientWaitSync(sync, 0, 0);
      if (status === gl.ALREADY_SIGNALED || status === gl.CONDITION_SATISFIED || status === gl.WAIT_FAILED) {
        gl.deleteSync(sync);
        resolve();
      } else {
        setTimeout(poll, 0);
      }
    };
    poll();
  });
}

function afterTwoFrames() {
  return new Promise(resolve => requestAnimationFrame(() => requestAnimationFrame(resolve)));
}
