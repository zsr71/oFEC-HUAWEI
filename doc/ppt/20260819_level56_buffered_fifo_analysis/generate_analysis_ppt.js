const pptxgen = require('pptxgenjs');

const pptx = new pptxgen();
pptx.layout = 'LAYOUT_WIDE';
pptx.author = 'oFEC';
pptx.company = 'oFEC';
pptx.subject = 'Eb/N0=3.1 dB, R_buf=32 的 Level 5/6 Buffered FIFO 运行分析';
pptx.title = 'Eb/N0=3.1 dB、R_buf=32：五六级 Buffered FIFO 运行分析';
pptx.lang = 'zh-CN';
pptx.theme = { headFontFace: 'Microsoft YaHei', bodyFontFace: 'Microsoft YaHei', lang: 'zh-CN' };

const C = {
  black: '111111', muted: '666666', gray: 'F2F2F2', gray2: 'E2E6EA', grid: 'B9C1C9',
  red: 'C00000', paleRed: 'FBECEC', blue: '476A86', paleBlue: 'EAF1F5',
  green: '4D7D5C', paleGreen: 'EDF5EF', orange: 'B66A17', paleOrange: 'FFF3E2',
  white: 'FFFFFF', dark: '27323A', darkBlue: '263D50', buffer: 'D9DEE3', l6: '5B7890', l5: '8DA6B8'
};
const FONT = 'Microsoft YaHei';

pptx.defineSlideMaster({
  title: 'MASTER',
  background: { color: C.white },
  objects: [
    { line: { x: 0.72, y: 0.96, w: 11.9, h: 0, line: { color: '252525', width: 1.3 } } },
    { text: { text: 'oFEC · Level 5/6 Buffered FIFO 运行分析', options: { x: 8.55, y: 7.12, w: 4.05, h: 0.18, fontFace: FONT, fontSize: 8, color: C.muted, align: 'right', margin: 0 } } },
  ],
  slideNumber: { x: 12.62, y: 7.12, color: C.muted, fontFace: FONT, fontSize: 8 },
});

function addText(s, text, x, y, w, h, opts = {}) {
  s.addText(text, { x, y, w, h, fontFace: FONT, fontSize: 14, color: C.black, margin: 0.04, fit: 'shrink', valign: 'mid', breakLine: false, ...opts });
}
function title(s, text, subtitle = '') {
  addText(s, text, 0.72, 0.40, 11.9, 0.4, { fontSize: 22, bold: true });
  if (subtitle) addText(s, subtitle, 0.72, 1.07, 11.7, 0.26, { fontSize: 12.5, color: C.muted });
}
function box(s, text, x, y, w, h, opts = {}) {
  s.addShape(pptx.ShapeType.rect, { x, y, w, h, rectRadius: 0.02, fill: { color: opts.fill || C.gray }, line: { color: opts.line || '777777', width: opts.lineWidth || 0.8 } });
  addText(s, text, x + 0.08, y + 0.05, w - 0.16, h - 0.1, { align: opts.align || 'center', valign: opts.valign || 'mid', bold: opts.bold || false, fontSize: opts.fontSize || 13, color: opts.color || C.black });
}
function keyBox(s, text, x, y, w, h, opts = {}) { box(s, text, x, y, w, h, { ...opts, line: opts.line || C.red, lineWidth: opts.lineWidth || 1.5, bold: true, fill: opts.fill || C.paleRed }); }
function arrow(s, x1, y1, x2, y2, color = '555555', width = 1.1) { s.addShape(pptx.ShapeType.line, { x: x1, y: y1, w: x2 - x1, h: y2 - y1, line: { color, width, endArrowType: 'triangle' } }); }
function bullets(s, items, x, y, w, h, fs = 14) {
  const runs = [];
  items.forEach((item) => runs.push({ text: item, options: { bullet: { indent: 12 }, hanging: 3, breakLine: true } }));
  addText(s, runs, x, y, w, h, { fontSize: fs, valign: 'top', paraSpaceAfterPt: 5 });
}
function table(s, rows, x, y, w, h, colW, fs = 12.2, rowH = 0.43) {
  s.addTable(rows, { x, y, w, h, border: { type: 'solid', color: C.grid, pt: 0.6 }, fill: C.white, color: C.black, fontFace: FONT, fontSize: fs, margin: 0.07, valign: 'mid', autoFit: false, rowH, colW, fillHeader: C.gray2, boldHeader: true });
}
function section(s, text, x, y, w, color = C.red) { addText(s, text, x, y, w, 0.25, { fontSize: 15, bold: true, color }); }
function metric(s, label, value, x, y, w, color = C.blue) {
  box(s, label, x, y, w, 0.38, { fill: C.gray2, line: C.grid, fontSize: 10.5, color: C.muted });
  addText(s, value, x, y + 0.40, w, 0.58, { fontSize: 21, bold: true, color, align: 'center' });
}
function memBar(s, x, y, w, h, start, label, fill, total = 76, color = C.black) {
  const segW = w * (22 / total);
  const sx = x + w * start / total;
  s.addShape(pptx.ShapeType.rect, { x: sx, y, w: segW, h, fill: { color: fill }, line: { color: C.white, width: 0.8 } });
  addText(s, label, sx + 0.03, y + 0.05, segW - 0.06, h - 0.1, { fontSize: 11.5, bold: true, color, align: 'center' });
}

// 1. Cover
{
  const s = pptx.addSlide('MASTER');
  s.background = { color: 'EEF1F3' };
  addText(s, 'oFEC 五 / 六级共享译码', 0.82, 1.15, 8.0, 0.35, { fontSize: 19, color: C.muted });
  addText(s, 'Eb/N0=3.1 dB、R_buf=32\n五六级 Buffered FIFO 运行分析', 0.82, 1.78, 8.7, 1.35, { fontSize: 31, bold: true, breakLine: true, valign: 'top' });
  addText(s, '从连续物理内存窗口运动，到 burst 下的队首延后与边界顶出', 0.86, 3.48, 9.5, 0.38, { fontSize: 17, color: C.blue });
  keyBox(s, 'Eb/N0 = 3.1 dB', 0.88, 4.50, 2.25, 0.62, { fontSize: 14 });
  box(s, 'R_buf = 32 行', 3.38, 4.50, 2.25, 0.62, { fontSize: 14 });
  box(s, '8 HISO + 8 SISO', 5.88, 4.50, 2.55, 0.62, { fontSize: 14 });
  addText(s, '基于 ofec_single 固定种子实测 · 2026-08-19', 0.88, 6.40, 5.8, 0.28, { fontSize: 12.5, color: C.muted });
  const bars = [0.55, 1.35, 1.98, 1.15, 0.52];
  bars.forEach((hh, i) => s.addShape(pptx.ShapeType.rect, { x: 9.20 + i * 0.62, y: 5.60 - hh, w: 0.39, h: hh, fill: { color: i === 2 ? C.red : C.blue }, line: { color: C.white, transparency: 100 } }));
  addText(s, 'burst', 9.25, 5.82, 2.6, 0.25, { fontSize: 11, color: C.muted, align: 'center' });
}

// 2. Run config and result
{
  const s = pptx.addSlide('MASTER'); title(s, '实验配置与总体结果', '本次展示只分析一个固定种子样本；结果用于解释运行机制，不替代多种子 BER 曲线');
  table(s, [
    [{ text: '配置项', options: { bold: true } }, { text: '本次设置', options: { bold: true } }, { text: '含义', options: { bold: true } }],
    ['Eb/N0', '3.1 dB', '固定单点运行'],
    ['R_buf', '32 block row', '可容纳 16 个“两行推进”时刻'],
    ['五六级共享资源', '8 HISO + 8 SISO', 'Level 5 / Level 6 共享服务池'],
    ['调度', 'Level5First + Group4LoadSortedMultiround', '沿用原五六级共享调度流程'],
    ['EarlyStop', 'condition=1, action=7', 'BCH + overall parity；沿用原 EarlyStopAction'],
    ['输入 / 输出节拍', '每个 t 固定两行', '窗口内部完成数不改变外部节拍'],
  ], 0.72, 1.48, 7.10, 3.45, [1.65, 2.45, 3.00], 11.4, 0.46);
  metric(s, 'Post-FEC BER', '2.09785×10⁻⁷', 8.35, 1.55, 1.85, C.red);
  metric(s, '统计错误数', '6 bit', 10.45, 1.55, 1.85, C.red);
  metric(s, '比较 bit 数', '28,600,704', 8.35, 3.05, 1.85, C.blue);
  metric(s, '服务时刻', '8,449', 10.45, 3.05, 1.85, C.blue);
  box(s, '观察重点\n不是“每个 t 是否完成”\n而是 burst 出现时\n窗口和 FIFO 如何共同吸收压力', 8.35, 4.65, 3.95, 1.35, { fill: C.paleBlue, line: C.blue, fontSize: 14, bold: true });
}

// 3. Physical memory layout
{
  const s = pptx.addSlide('MASTER'); title(s, '连续物理内存：R_buf=32 时的 76 行布局');
  addText(s, '这里的行是 block row；列方向的 8 个 block 不在本页展开。窗口运动只改变起始行 S_t，Level 5 / Level 6 相对位置始终保持。', 0.75, 1.30, 11.3, 0.42, { fontSize: 14, color: C.muted });
  const x = 0.95, y = 2.18, w = 10.85, h = 0.92;
  const segs = [{ n: 32, label: 'Buffer\n[0,32)', fill: C.buffer, color: C.black }, { n: 22, label: 'Level 6\n[32,54)', fill: C.l6, color: C.white }, { n: 22, label: 'Level 5\n[54,76)', fill: C.l5, color: C.black }];
  let cur = x; segs.forEach((a) => { const sw = w * a.n / 76; s.addShape(pptx.ShapeType.rect, { x: cur, y, w: sw, h, fill: { color: a.fill }, line: { color: C.white, width: 1 } }); addText(s, a.label, cur, y + 0.08, sw, 0.68, { fontSize: 15, bold: true, color: a.color, align: 'center' }); cur += sw; });
  addText(s, '基准窗口：S_t = 32', 0.98, 3.55, 3.2, 0.3, { fontSize: 16, bold: true, color: C.blue });
  addText(s, '当前五六级解码区覆盖 [32,76)，共 44 行', 4.05, 3.55, 5.4, 0.3, { fontSize: 15, color: C.muted });
  section(s, '窗口起点变化示例', 0.95, 4.28, 2.4);
  const states = [{ s: 32, text: 'S=32\n基准位置' }, { s: 24, text: 'S=24\n退后 8 行' }, { s: 0, text: 'S=0\nbuffer 用尽' }];
  states.forEach((st, i) => { const yy = 4.78 + i * 0.56; addText(s, st.text, 1.0 + i * 3.92, yy, 1.45, 0.4, { fontSize: 13, bold: true, color: i === 2 ? C.red : C.black, align: 'center' }); const bx = 2.55 + i * 3.92; const by = yy + 0.03; const bw = 2.45; s.addShape(pptx.ShapeType.rect, { x: bx, y: by, w: bw, h: 0.31, fill: { color: C.buffer }, line: { color: C.grid, width: 0.5 } }); s.addShape(pptx.ShapeType.rect, { x: bx + bw * st.s / 76, y: by, w: bw * 22 / 76, h: 0.31, fill: { color: C.l6 }, line: { color: C.white, width: 0.5 } }); s.addShape(pptx.ShapeType.rect, { x: bx + bw * (st.s + 22) / 76, y: by, w: bw * 22 / 76, h: 0.31, fill: { color: C.l5 }, line: { color: C.white, width: 0.5 } }); if (i < 2) arrow(s, bx + bw + 0.15, by + 0.15, bx + bw + 0.65, by + 0.15, C.blue); });
  addText(s, '注意：S_t 表示窗口的逻辑/物理位置；它不是 FIFO 深度，也不是输出行数。', 0.95, 6.58, 11.0, 0.3, { fontSize: 14, color: C.red, bold: true });
}

// 4. Per-time state machine
{
  const s = pptx.addSlide('MASTER'); title(s, '每个 t 的固定节拍与内部服务顺序');
  const steps = [['①', '底部进入两行', '新 batch B_t 到达'], ['②', '队首服务', '只服务 FIFO 队首；队首未完成时不越过'], ['③', '顶部输出两行', '物理 SRAM 节拍固定'], ['④', '更新 S_t', '由本时刻内部完成数 C_t 决定']];
  steps.forEach((st, i) => { const x = 0.78 + i * 3.08; keyBox(s, st[0], x, 1.55, 0.45, 0.45, { fontSize: 14, fill: C.white }); box(s, st[1], x + 0.62, 1.48, 2.12, 0.62, { fill: i === 1 ? C.paleRed : C.gray, line: i === 1 ? C.red : C.grid, bold: true, fontSize: 14 }); addText(s, st[2], x + 0.62, 2.25, 2.12, 0.35, { fontSize: 12, color: C.muted, align: 'center' }); if (i < steps.length - 1) arrow(s, x + 2.83, 1.79, x + 3.02, 1.79); });
  keyBox(s, '物理层不变量：每个 t 始终输出 2 个 block row', 2.05, 3.15, 8.8, 0.72, { fontSize: 18, fill: C.paleBlue, line: C.blue });
  const cases = [
    ['C_t=0', 'S: -2', '队首 Bx 不完成，下一时刻仍是 Bx', C.red, C.paleRed],
    ['C_t=1', 'S: 不变', '队首 Bx 完成，下一时刻切换 Bx+1', C.blue, C.paleBlue],
    ['C_t≥2', 'S: 恢复', '普通完成 + 全 EarlyStop 快路径，消化积压', C.green, C.paleGreen],
  ];
  cases.forEach((c, i) => { const x = 0.82 + i * 4.05; box(s, c[0], x, 4.42, 1.18, 0.58, { fill: c[4], line: c[3], bold: true, fontSize: 14 }); box(s, c[1], x + 1.35, 4.42, 1.15, 0.58, { fill: C.white, line: c[3], fontSize: 14 }); addText(s, c[2], x, 5.28, 3.55, 0.7, { fontSize: 13, color: C.muted, valign: 'top' }); });
  addText(s, '这里的 C_t 是五六级内部完成的 batch 数，不是对外一次输出了几个 batch。', 0.9, 6.35, 11.2, 0.34, { fontSize: 14, bold: true, color: C.red });
}

// 5. FIFO depth and queue head
{
  const s = pptx.addSlide('MASTER'); title(s, 'FIFO depth 与队首：两个容易混淆的量');
  keyBox(s, 'FIFO depth = 当前 FIFO 中保存的 64-code batch 数量', 1.55, 1.48, 9.55, 0.72, { fontSize: 18 });
  table(s, [
    [{ text: '字段', options: { bold: true } }, { text: '本页定义', options: { bold: true } }, { text: '不是', options: { bold: true } }],
    ['fifo_depth_before', '新 batch 到达后、服务开始前的 batch 数', '不是 block row 数'],
    ['fifo_depth_after', '本时刻退休/顶出后剩余的 batch 数', '不是本时刻输出行数'],
    ['head_batch', 'FIFO 当前最老、下一次优先服务的 batch', '不会因 C_t=0 自动变化'],
    ['forced_evicted_global_rows', '被顶出 batch 中仍未完成 code 的诊断行列表', '不是物理顶出行数'],
  ], 0.85, 2.62, 11.35, 2.35, [2.55, 5.25, 3.55], 12.3, 0.48);
  const qx = 1.0, qy = 5.42; ['B2985', 'B2986', 'B2987', '…', 'B3001'].forEach((b, i) => { box(s, b, qx + i * 1.58, qy, 1.22, 0.54, { fill: i === 0 ? C.paleRed : C.paleBlue, line: i === 0 ? C.red : C.blue, bold: i === 0, fontSize: 13 }); if (i < 4) arrow(s, qx + 1.28 + i * 1.58, qy + 0.27, qx + 1.48 + i * 1.58, qy + 0.27); });
  addText(s, '队首是 B2985 时，如果 C_t=0，下一时刻仍然服务 B2985；只有 B2985 完成或被顶出后，B2986 才成为队首。', 1.0, 6.22, 11.0, 0.42, { fontSize: 14, color: C.red, bold: true });
}

// 6. Whole run stats
{
  const s = pptx.addSlide('MASTER'); title(s, '全帧运行：大多数时刻稳定，burst 时出现积压');
  const stats = [['总服务时刻', '8,449', C.blue], ['C_t=0', '61', C.red], ['C_t=1', '8,346', C.blue], ['C_t≥2', '42', C.green], ['ForcedEvicted batch', '17', C.red], ['全 EarlyStop batch', '189', C.green]];
  stats.forEach((a, i) => metric(s, a[0], a[1], 0.85 + (i % 3) * 4.05, 1.50 + Math.floor(i / 3) * 1.30, 3.35, a[2]));
  section(s, 'S_t 的覆盖情况', 0.9, 4.28, 2.4);
  const dist = [{ s: 0, n: 234, fill: C.red }, { s: 2, n: 82, fill: 'D8815A' }, { s: 4, n: 167, fill: 'D8A16E' }, { s: 6, n: 204, fill: 'D6BC86' }, { s: 8, n: 50, fill: 'C8C98D' }, { s: 10, n: 5, fill: 'A7C38A' }, { s: 12, n: 8, fill: '8BB184' }, { s: 14, n: 48, fill: '72A27E' }, { s: 16, n: 11, fill: '60967B' }, { s: 18, n: 16, fill: '5D8C78' }, { s: 20, n: 79, fill: '557F76' }, { s: 22, n: 197, fill: '52747A' }, { s: 24, n: 278, fill: '4C6A78' }, { s: 26, n: 232, fill: '4A6075' }, { s: 28, n: 319, fill: '485772' }, { s: 30, n: 780, fill: '44506D' }, { s: 32, n: 5739, fill: C.blue }];
  const max = 5739, chartX = 1.0, chartY = 4.88, chartW = 10.8, chartH = 1.15; dist.forEach((d, i) => { const bw = chartW / dist.length * 0.82; const xx = chartX + i * chartW / dist.length + 0.06; const hh = chartH * d.n / max; s.addShape(pptx.ShapeType.rect, { x: xx, y: chartY + chartH - hh, w: bw, h: hh, fill: { color: d.fill }, line: { color: C.white, transparency: 100 } }); addText(s, String(d.s), xx - 0.05, chartY + chartH + 0.05, bw + 0.10, 0.22, { fontSize: 8.5, color: C.muted, align: 'center' }); });
  addText(s, 'S_t', 0.72, 5.90, 0.25, 0.2, { fontSize: 10, color: C.muted });
  addText(s, '基准位置 S=32 占 5,739 / 8,449 = 67.9%；buffer 只在 burst 时段被逐步占用。', 1.0, 6.32, 11.0, 0.34, { fontSize: 14, color: C.blue, bold: true });
}

// 7. Correct burst timeline
{
  const s = pptx.addSlide('MASTER'); title(s, '真实 burst 片段：S=26 退到 S=0 的完整因果链', '下面只抽取关键时刻；省略的 C_t=1 时刻会完成当前队首并切换到下一个 batch');
  table(s, [
    [{ text: 't', options: { bold: true } }, { text: '队首', options: { bold: true } }, { text: 'C_t', options: { bold: true } }, { text: 'S_t→S_{t+1}', options: { bold: true } }, { text: '真实含义', options: { bold: true } }],
    ['2897', 'B2894', '0', '26→24', 'B2894 未完成；下一时刻仍是 B2894'],
    ['2898', 'B2894', '1', '24→24', 'B2894 完成；之后才换 B2895'],
    ['2908', 'B2904', '0', '24→22', 'B2904 未完成；下一时刻仍是 B2904'],
    ['2909', 'B2904', '1', '22→22', 'B2904 完成；之后才换 B2905'],
    ['2910', 'B2905', '0', '22→20', 'B2905 未完成；下一时刻仍是 B2905'],
    ['2911', 'B2905', '1', '20→20', 'B2905 完成；之后才换 B2906'],
    ['2999', 'B2984', '0', '2→0', '窗口到达边界，B2984 仍留在队首'],
    ['3000', 'B2984', '1', '0→0', 'B2984 完成，但窗口不能再恢复'],
    ['3001', 'B2985', '0', '0→0', 'B2985 服务后仍未完成，ForcedEvicted'],
  ], 0.65, 1.43, 12.05, 4.25, [0.72, 1.20, 0.72, 1.55, 7.86], 11.1, 0.43);
  keyBox(s, '关键规则：C_t=0 时，S_t 后退，但队首不变；队首只有完成或被顶出后才切换。', 1.05, 5.95, 10.9, 0.70, { fontSize: 16, fill: C.paleRed });
}

// 8. S=0 forced eviction
{
  const s = pptx.addSlide('MASTER'); title(s, 'S_t=0 时的 ForcedEvicted：保持两行节拍的边界保护');
  addText(s, '本次共 17 个 batch 触发 ForcedEvicted。它们都在 S_t=0 时先获得一次普通服务机会，仍未完成后才被强制退休。', 0.78, 1.28, 11.4, 0.38, { fontSize: 14, color: C.muted });
  table(s, [
    [{ text: 't', options: { bold: true } }, { text: '队首 batch', options: { bold: true } }, { text: 'FIFO depth', options: { bold: true } }, { text: '物理输出', options: { bold: true } }, { text: '诊断信息', options: { bold: true } }],
    ['3001', 'B2985', '17→16 batches', '固定 2 block row', '2 个 code 仍未完成：96219, 96223'],
    ['3003', 'B2987', '17→16 batches', '固定 2 block row', '3 个 code 仍未完成：96277, 96287, 95917'],
    ['3005', 'B2989', '17→16 batches', '固定 2 block row', '1 个 code 仍未完成：96351'],
    ['3006', 'B2990', '17→16 batches', '固定 2 block row', '2 个 code 仍未完成：96379, 96018'],
    ['3549', 'B3533', '17→16 batches', '固定 2 block row', '1 个 code 仍未完成：113743'],
    ['3551', 'B3535', '17→16 batches', '固定 2 block row', '1 个 code 仍未完成：113810'],
  ], 0.65, 1.83, 12.05, 2.92, [0.75, 1.32, 1.58, 1.68, 6.72], 11.2, 0.46);
  const qx = 1.05, qy = 5.20;
  box(s, 'S=0', qx, qy, 0.95, 0.55, { fill: C.paleRed, line: C.red, bold: true, fontSize: 15 });
  arrow(s, qx + 1.15, qy + 0.28, qx + 2.0, qy + 0.28, C.red);
  box(s, '普通服务仍未完成', qx + 2.10, qy, 2.25, 0.55, { fill: C.paleRed, line: C.red, fontSize: 13 });
  arrow(s, qx + 4.55, qy + 0.28, qx + 5.40, qy + 0.28, C.red);
  keyBox(s, 'ForcedEvicted\n退休一个队首 batch', qx + 5.50, qy - 0.10, 2.45, 0.75, { fontSize: 13 });
  arrow(s, qx + 8.15, qy + 0.28, qx + 9.0, qy + 0.28, C.blue);
  box(s, '下一时刻继续输出两行\nFIFO 维持可服务节拍', qx + 9.10, qy - 0.02, 2.45, 0.62, { fill: C.paleBlue, line: C.blue, fontSize: 12.5 });
  addText(s, '物理输出两行/t 与 ForcedEvicted batch 是两个不同概念；forced_evicted_global_rows 只是该 batch 未完成 code 的行列表。', 0.95, 6.32, 11.2, 0.34, { fontSize: 13.5, color: C.red, bold: true });
}

// 9. Physical SRAM snapshots through a burst
{
  const s = pptx.addSlide('MASTER'); title(s, 'burst 的三个物理内存状态：基准 → 回退 → 饱和', '每幅图都是同一块 76-row 连续 SRAM；这里只改变解码窗口位置 S_t，不改变 Level 6 / Level 5 的相对间距');
  const cards = [
    { x: 0.60, stage: '阶段 1：基准运行', st: 'S_t=32', head: '队首：B_x 当期完成', fifo: 'FIFO：约 0–1 batch', note: 'C_t=1\n输入推进与完成抵消', start: 32, accent: C.blue, pale: C.paleBlue },
    { x: 4.48, stage: '阶段 2：burst 积压', st: 'S_t=16', head: '队首：B_x pending', fifo: 'FIFO：持续累积', note: 'C_t=0\n队首保持，窗口退 2 行', start: 16, accent: C.orange, pale: C.paleOrange },
    { x: 8.36, stage: '阶段 3：buffer 用尽', st: 'S_t=0', head: '队首：B_x 仍 pending', fifo: 'FIFO：满载约 16 batch', note: 'C_t=0\n不能再退，必要时 ForcedEvicted', start: 0, accent: C.red, pale: C.paleRed },
  ];
  cards.forEach((a, i) => {
    const x = a.x, y = 1.56, w = 3.45;
    keyBox(s, a.stage, x, y, w, 0.48, { fontSize: 13.5, fill: a.pale, line: a.accent });
    addText(s, a.st, x, y + 0.60, w, 0.30, { fontSize: 18, bold: true, color: a.accent, align: 'center' });
    // A vertical physical-memory view: top is physical output side, bottom is new-input side.
    const mx = x + 0.45, my = 2.24, mw = 2.55, mh = 2.42;
    s.addShape(pptx.ShapeType.rect, { x: mx, y: my, w: mw, h: mh, fill: { color: C.white }, line: { color: C.dark, width: 1.1 } });
    const unit = mh / 76;
    const l6y = my + a.start * unit;
    const l5y = my + (a.start + 22) * unit;
    s.addShape(pptx.ShapeType.rect, { x: mx, y: l6y, w: mw, h: 22 * unit, fill: { color: C.l6 }, line: { color: C.white, width: 0.6 } });
    s.addShape(pptx.ShapeType.rect, { x: mx, y: l5y, w: mw, h: 22 * unit, fill: { color: C.l5 }, line: { color: C.white, width: 0.6 } });
    addText(s, 'L6\n[' + a.start + ',' + (a.start + 22) + ')', mx + 0.04, l6y + 0.01, mw - 0.08, Math.max(0.26, 22 * unit - 0.02), { fontSize: 11, bold: true, color: C.white, align: 'center' });
    addText(s, 'L5\n[' + (a.start + 22) + ',' + (a.start + 44) + ')', mx + 0.04, l5y + 0.01, mw - 0.08, Math.max(0.26, 22 * unit - 0.02), { fontSize: 11, bold: true, color: C.black, align: 'center' });
    if (a.start > 0) {
      addText(s, '可用 buffer\n' + a.start + ' 行', mx + 0.08, my + 0.05, mw - 0.16, Math.max(0.22, a.start * unit - 0.05), { fontSize: 10.5, color: C.muted, align: 'center' });
    } else {
      addText(s, 'buffer\n已用尽', mx + 0.10, my + 0.04, mw - 0.20, 0.25, { fontSize: 10.5, color: C.red, bold: true, align: 'center' });
    }
    addText(s, '顶部：固定输出 2 行', x + 0.12, my - 0.26, 3.20, 0.20, { fontSize: 10.5, color: C.muted, align: 'center' });
    addText(s, '底部：固定进入 2 行', x + 0.12, my + mh + 0.06, 3.20, 0.20, { fontSize: 10.5, color: C.muted, align: 'center' });
    box(s, a.head, x + 0.08, 5.12, w - 0.16, 0.43, { fill: C.white, line: a.accent, fontSize: 12.2, bold: true });
    box(s, a.fifo, x + 0.08, 5.64, w - 0.16, 0.40, { fill: C.gray, line: C.grid, fontSize: 11.8 });
    box(s, a.note, x + 0.08, 6.14, w - 0.16, 0.56, { fill: a.pale, line: a.accent, fontSize: 11.5, color: a.accent, bold: true });
    if (i < cards.length - 1) {
      arrow(s, x + w + 0.03, 3.46, x + w + 0.35, 3.46, a.accent, 1.5);
    }
  });
  addText(s, '读图顺序：burst 中多次 C_t=0 使 S_t 逐步从 32 减少到 0；每次 C_t=0 时，当前队首 batch 仍保持不变，下一时刻继续获得服务。', 0.75, 6.88, 11.7, 0.26, { fontSize: 12.7, color: C.black, bold: true, align: 'center' });
}

// 10. R_buf=32 vs 64 summary
{
  const s = pptx.addSlide('MASTER'); title(s, 'R_buf=32 vs 64：缓冲加倍后，边界顶出消失', '两组使用相同 Eb/N0=3.1 dB、种子、共享资源与调度；唯一改变的是 Buffer 长度');
  table(s, [
    [{ text: '指标', options: { bold: true } }, { text: 'R_buf=32', options: { bold: true } }, { text: 'R_buf=64', options: { bold: true } }, { text: '解释', options: { bold: true } }],
    ['总 SRAM 高度', '76 行 = 32+22+22', '108 行 = 64+22+22', 'Level 5/6 解码区不变，只扩展 buffer'],
    ['最大额外等待空间', '16 t', '32 t', '每 t 固定推进 2 block row'],
    ['最小 S_t', '0', '6', '64 行下未触及物理边界'],
    ['ForcedEvicted batch', '17', '0', '原先被截断的 batch 均可继续等待'],
    ['最大 FIFO depth（服务前）', '17 batch', '30 batch', '更大 buffer 吸收更深积压'],
    ['正常 batch 平均等待', '1.885 t', '4.406 t', '延迟增加，换取完成率'],
    ['最长正常等待', '16 t', '29 t', '仍未达到 R_buf=64 的 32 t 上限'],
    ['Post-FEC BER', '6 / 28,600,704\n2.09785×10⁻⁷', '5 / 28,600,704\n1.74821×10⁻⁷', '单种子下少 1 个错误；方向正确但统计量很小'],
  ], 0.55, 1.38, 12.30, 4.55, [2.10, 2.15, 2.15, 5.90], 11.1, 0.47);
  keyBox(s, '这次 burst 的最大实际等待为 29 t，因此 R_buf=32 不够，而 R_buf=64 仍保留 3 t / 6 行余量。', 1.00, 6.15, 11.10, 0.60, { fontSize: 15.5, fill: C.paleBlue, line: C.blue });
}

// 11. R_buf=32 vs 64 case comparison
{
  const s = pptx.addSlide('MASTER'); title(s, '同一批次在 32 / 64 行 buffer 下的不同命运', 'B2985、B3535 是最有代表性的两个案例：R_buf=32 触边顶出；R_buf=64 继续等待后 Normal 完成');
  const cols = [
    { x: 0.70, name: 'B2985', arrive: '到达 t=2985', r32: ['t=2999：B2984 使 S:2→0', 't=3001：B2985 服务后仍 pending', '等待 16 t，ForcedEvicted'], r64: ['t=3001：B2985 pending，S:32→30', 't=3002：B2985 Normal 完成', '等待 17 t，未触边'] },
    { x: 6.80, name: 'B3535', arrive: '到达 t=3535', r32: ['t=3549：B3533 ForcedEvicted', 't=3551：B3535 服务后仍 pending', '等待 16 t，ForcedEvicted'], r64: ['t=3563：B3535 pending，S:8→6', 't=3564：B3535 Normal 完成', '等待 29 t，仍有 6 行余量'] },
  ];
  cols.forEach((a) => {
    keyBox(s, a.name, a.x, 1.46, 5.15, 0.56, { fontSize: 17, fill: C.white });
    addText(s, a.arrive, a.x, 2.14, 5.15, 0.24, { fontSize: 13.5, color: C.muted, align: 'center' });
    box(s, 'R_buf=32', a.x, 2.70, 2.35, 0.48, { fill: C.paleRed, line: C.red, bold: true, fontSize: 14, color: C.red });
    box(s, 'R_buf=64', a.x + 2.80, 2.70, 2.35, 0.48, { fill: C.paleGreen, line: C.green, bold: true, fontSize: 14, color: C.green });
    bullets(s, a.r32, a.x + 0.04, 3.36, 2.26, 1.45, 12.5);
    bullets(s, a.r64, a.x + 2.84, 3.36, 2.26, 1.45, 12.5);
    box(s, '结局：强制退休\n不能再获得后续服务', a.x, 5.05, 2.35, 0.76, { fill: C.paleRed, line: C.red, fontSize: 13, bold: true, color: C.red });
    box(s, '结局：Normal 完成\n17 个原顶出 batch 全部恢复', a.x + 2.80, 5.05, 2.35, 0.76, { fill: C.paleGreen, line: C.green, fontSize: 13, bold: true, color: C.green });
  });
  addText(s, '机制结论：R_buf=64 没有改变 FIFO 顺序或单个 t 的资源预算；它只额外提供了继续服务同一队首 batch 的时间，避免 S=0 边界截断。', 0.80, 6.28, 11.75, 0.38, { fontSize: 14.5, bold: true, color: C.blue, fill: { color: C.paleBlue }, margin: 0.10, align: 'center' });
}

// 12. What is measured / implementation boundary
{
  const s = pptx.addSlide('MASTER'); title(s, '本次结果能说明什么？当前软件仿真的边界在哪里？');
  const cols = [
    ['本次可以直接确认', ['FIFO 队首保持与完成后切换关系', 'C_t=0 / 1 / ≥2 对 S_t 的影响', 'R_buf=32 对最大等待时刻的约束', 'S=0 时 ForcedEvicted 的触发条件']],
    ['不能直接等同于', ['真实 SRAM 每个 block row 已经在软件中整体搬移', '所有顶出行都对应某一个 ForcedEvicted code', '单种子 BER 就代表总体性能', 'ForcedEvicted 次数等于物理顶出次数']],
    ['后续实现核对点', ['把连续 SRAM 的地址映射接入 batch 读写', '补充 FIFO drain 或明确帧尾统计策略', '对 R_buf 扫描统一统计口径', '区分物理两行输出与 batch 退休统计']],
  ];
  cols.forEach((c, i) => { const x = 0.76 + i * 4.08; keyBox(s, c[0], x, 1.55, 3.45, 0.62, { fontSize: 14, fill: i === 0 ? C.paleBlue : i === 1 ? C.paleRed : C.paleOrange, line: i === 0 ? C.blue : i === 1 ? C.red : C.orange }); bullets(s, c[1], x + 0.12, 2.45, 3.18, 2.5, 13.2); });
  addText(s, '这页的核心态度：调度状态序列已经能解释 burst 的窗口行为；物理 SRAM 的真实搬移仍需要单独做代码级验证。', 0.9, 5.75, 11.1, 0.62, { fontSize: 16, bold: true, color: C.blue, fill: { color: C.paleBlue }, margin: 0.13 });
}

// 13. Takeaways
{
  const s = pptx.addSlide('MASTER'); title(s, '结论：Buffered FIFO 吸收短时 burst，但不能突破长期服务能力');
  const take = [
    ['1', '大多数时刻稳定', '本次 8,449 个服务时刻中，8,346 个时刻 C_t=1，窗口保持节拍稳定。'],
    ['2', 'burst 时先退窗口', 'C_t=0 时队首不变，S_t 每次退 2 行，FIFO 继续按先进先出等待。'],
    ['3', '低负载时再追回', 'C_t≥2 时使用全 EarlyStop 快路径逐步恢复 S_t，窗口可回到基准位置。'],
    ['4', '边界时优先保吞吐', 'S_t=0 且队首仍未完成时，ForcedEvicted 保证物理输出仍是 2 行/t。'],
  ];
  take.forEach((a, i) => { const y = 1.48 + i * 0.91; keyBox(s, a[0], 0.95, y, 0.48, 0.48, { fontSize: 14, fill: C.white }); addText(s, a[1], 1.70, y + 0.02, 2.4, 0.30, { fontSize: 16, bold: true, color: i === 3 ? C.red : C.blue }); addText(s, a[2], 4.10, y + 0.02, 7.8, 0.40, { fontSize: 14, color: C.black }); });
  keyBox(s, '一句话结论', 0.95, 5.42, 1.65, 0.62, { fontSize: 15, fill: C.paleRed });
  addText(s, 'Buffered FIFO 把“某个时刻处理不完”转换为“后续时刻继续处理”；它解决短时突发，不改变长期平均吞吐上限。', 2.85, 5.48, 9.3, 0.55, { fontSize: 18, bold: true, color: C.black });
  addText(s, '展示建议：重点讲第 4、5、7、8、10、11 页；第 12 页用于主动说明当前仿真边界。', 0.95, 6.45, 11.0, 0.28, { fontSize: 12.5, color: C.muted });
}

pptx.writeFile({ fileName: __dirname + '/level56_buffered_fifo_ebn03p1_rbuf32_analysis.pptx' });
