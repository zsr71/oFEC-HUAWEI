const pptxgen = require('pptxgenjs');

const pptx = new pptxgen();
pptx.layout = 'LAYOUT_WIDE';
pptx.author = 'oFEC';
pptx.subject = 'Editable PowerPoint sample';
pptx.title = 'oFEC Level 5/6 editable sample';
pptx.company = 'oFEC';
pptx.lang = 'zh-CN';
pptx.theme = {
  headFontFace: 'Microsoft YaHei',
  bodyFontFace: 'Microsoft YaHei',
  lang: 'zh-CN',
};
pptx.defineSlideMaster({
  title: 'MASTER',
  background: { color: 'FFFFFF' },
  objects: [
    { line: { x: 0.75, y: 0.96, w: 11.85, h: 0, line: { color: '222222', width: 1.5 } } },
    { text: { text: 'oFEC · Level 5/6 早停更新与 History 修正', options: { x: 9.55, y: 7.12, w: 3.0, h: 0.18, fontFace: 'Microsoft YaHei', fontSize: 8, color: '666666', align: 'right', margin: 0 } } },
  ],
  slideNumber: { x: 12.62, y: 7.12, color: '666666', fontFace: 'Microsoft YaHei', fontSize: 8 },
});

const C = { black: '111111', gray: 'F2F2F2', grid: 'BFC5CC', red: 'C00000', muted: '666666', white: 'FFFFFF' };
const FONT = 'Microsoft YaHei';

function title(slide, text) {
  slide.addText(text, { x: 0.75, y: 0.42, w: 11.9, h: 0.38, fontFace: FONT, fontSize: 22, bold: true, color: C.black, margin: 0 });
}

function addText(slide, text, x, y, w, h, opts = {}) {
  slide.addText(text, { x, y, w, h, fontFace: FONT, fontSize: 14, color: C.black, margin: 0.04, breakLine: false, fit: 'shrink', valign: 'mid', ...opts });
}

function box(slide, text, x, y, w, h, opts = {}) {
  slide.addShape(pptx.ShapeType.rect, { x, y, w, h, rectRadius: 0, fill: { color: opts.fill || C.gray }, line: { color: opts.line || '777777', width: opts.lineWidth || 0.8 } });
  addText(slide, text, x + 0.08, y + 0.06, w - 0.16, h - 0.12, { align: 'center', valign: 'mid', bold: opts.bold || false, fontSize: opts.fontSize || 13, color: opts.color || C.black });
}

function keyBox(slide, text, x, y, w, h) {
  box(slide, text, x, y, w, h, { line: C.red, lineWidth: 1.5, bold: true });
}

function arrow(slide, x1, y1, x2, y2) {
  slide.addShape(pptx.ShapeType.line, { x: x1, y: y1, w: x2 - x1, h: y2 - y1, line: { color: '555555', width: 1, endArrowType: 'triangle' } });
}

function bullets(slide, items, x, y, w, h, fontSize = 14) {
  addText(slide, items.map((item) => ({ text: item, options: { bullet: { indent: 12 }, hanging: 3, breakLine: true } })), x, y, w, h, { fontSize, valign: 'top', paraSpaceAfterPt: 5 });
}

// Page 3: flow information represented entirely by editable table cells.
{
  const s = pptx.addSlide('MASTER');
  title(s, '三种方案的总体差异');
  addText(s, '公共流程和模式差异合并如下：', 0.75, 1.12, 5.3, 0.3, { fontSize: 15 });
  s.addTable([
    [{ text: '流程阶段', options: { bold: true } }, { text: '统一执行内容', options: { bold: true } }, { text: '三种模式的差异', options: { bold: true } }],
    ['1. 行级预处理', '检测 EarlyStop，并完成 Hybrid 分类', '三种模式相同'],
    ['2. 负载统计', '统计 16 个 group 的普通候选，计算非空组数 K', '三种模式相同'],
    ['3. 普通调度', '按负载排序，分轮使用最多 8 次 group entry', '三种模式相同'],
    ['4. 进入记录', '记录每组是否普通调度进入：group_entered', '三种模式相同'],
    ['5. EarlyStop 更新', '命中行决定是否写入新的外信息', 'AllGroups：所有命中组；EnteredGroupsOnly：仅普通进入组；FillIdleEntries：普通进入组加空闲补位组'],
    ['6. History 写回', '有新外信息则合并；无新外信息则透传当前 prior', '三种模式共用修正后的 History 语义'],
  ], { x: 0.75, y: 1.52, w: 11.85, h: 4.95, border: { type: 'solid', color: C.grid, pt: 0.6 }, fill: C.white, color: C.black, fontFace: FONT, fontSize: 12.2, margin: 0.07, valign: 'mid', autoFit: false, rowH: 0.64, colW: [1.58, 3.35, 6.92], bold: false, fillHeader: 'E3E3E3', boldHeader: true });
  addText(s, '三种模式只改变 EarlyStop 外信息的更新资格，普通 HISO/SISO 的 core 和 MUX 规则不变。', 0.75, 6.63, 11.8, 0.3, { fontSize: 13.5 });
}

// Page 8: explicit distinction between resource state and strategy.
{
  const s = pptx.addSlide('MASTER');
  title(s, '空闲 entry 的具体利用方式');
  box(s, '第一轮\n4 组', 1.05, 1.55, 1.45, 0.75, { fontSize: 15 });
  box(s, '第二轮\n2 组', 3.0, 1.55, 1.45, 0.75, { fontSize: 15 });
  box(s, 'used = 6\nremaining = 2', 4.95, 1.55, 2.05, 0.75, { fontSize: 15 });
  arrow(s, 2.5, 1.93, 3.0, 1.93); arrow(s, 4.45, 1.93, 4.95, 1.93);
  keyBox(s, 'FillIdleEntries\n补 2 个未进入组', 8.35, 2.35, 2.4, 0.78); box(s, 'EnteredGroupsOnly\n2 个 entry 浪费', 8.38, 1.15, 2.34, 0.75, { fontSize: 13 });
  arrow(s, 7.0, 1.72, 8.35, 2.7); arrow(s, 7.0, 1.72, 8.38, 1.52);
  addText(s, '示例：第一轮进入 4 个 group，第二轮进入 2 个 group。\n这里“状态”只描述资源使用情况，“方案”描述三种策略如何处理一个状态。', 0.75, 3.45, 11.8, 0.65, { fontSize: 14, valign: 'top' });
  s.addTable([
    [{ text: '类别', options: { bold: true } }, { text: '名称', options: { bold: true } }, { text: '动作', options: { bold: true } }],
    ['调度状态', 'used=6，remaining=2', '普通调度已使用 6 次，还剩 2 次 entry'],
    ['方案', 'EnteredGroupsOnly', '剩余 2 个 entry 不使用，只更新已进入组'],
    ['方案', 'FillIdleEntries', '用 2 个空闲 entry 补两个尚未进入的最小 index group'],
    ['方案', 'AllGroups', '不受 group 是否进入限制，所有 EarlyStop 命中组都可更新'],
  ], { x: 2.55, y: 4.35, w: 8.3, h: 2.1, border: { type: 'solid', color: C.grid, pt: 0.6 }, fontFace: FONT, fontSize: 12.4, margin: 0.07, valign: 'mid', rowH: 0.42, colW: [1.25, 2.25, 4.8], fillHeader: 'E3E3E3', boldHeader: true });
  addText(s, 'FillIdleEntries 恢复的是部分 EarlyStop 更新机会，不增加 HISO/SISO 解码容量。', 0.75, 6.65, 11.8, 0.3, { fontSize: 13.5 });
}

// Page 11: the root cause as editable text, without a flowchart.
{
  const s = pptx.addSlide('MASTER');
  title(s, '根因：produced=false 时旧 History 被保留');
  addText(s, '这里 produced 表示“当前行是否真正产生了新的译码输出”，不是“这一行是否存在”。', 0.75, 1.2, 11.7, 0.35, { fontSize: 15 });
  addText(s, '修复前：', 0.75, 1.85, 1.4, 0.3, { fontSize: 15, bold: true });
  bullets(s, ['produced=true：history = prior + new_extrinsic', 'produced=false：不写 history，保留旧窗口值或初始化 0'], 0.95, 2.2, 10.8, 0.85, 15);
  box(s, '错误输出\nfinal_output_old = channel + stale_or_zero_history', 1.0, 3.35, 5.2, 0.95, { line: C.red, lineWidth: 1.4, bold: true, fontSize: 14 });
  box(s, '正确应为\nfinal_output = channel + 当前 prior', 7.05, 3.35, 4.9, 0.95, { line: '3F7F5F', lineWidth: 1.2, bold: true, fontSize: 14 });
  addText(s, '因此，produced=false 的含义是“本行不新增信息”，而不是“本行的已有信息清零”。', 0.75, 4.75, 11.8, 0.35, { fontSize: 15 });
  addText(s, '受影响范围包括资源不足的 Unscheduled 行，以及未获得 EarlyStop 更新资格的命中行。', 0.75, 5.2, 11.8, 0.35, { fontSize: 15 });
  addText(s, '根因结论：不是 EarlyStop 判决本身错误，而是 produced=false 时没有把当前 prior 透传到最后一个 tile 的 History。', 0.95, 6.05, 11.2, 0.55, { fontSize: 15, color: '536B82', fill: { color: 'F3F3F3' }, margin: 0.12, breakLine: false });
}

// Page 12: editable rectangular branch flow.
{
  const s = pptx.addSlide('MASTER');
  title(s, 'History 修正：透传当前 Prior');
  box(s, '当前 tile 输入\n已有 prior', 0.85, 1.8, 1.7, 0.85, { fontSize: 14 });
  keyBox(s, '是否产生新的\n外信息？', 3.0, 1.8, 2.0, 0.85);
  box(s, '合并当前结果\nhistory = prior + new_extrinsic', 6.0, 1.15, 2.8, 1.0, { fontSize: 13 });
  keyBox(s, '不新增译码结果\n直接透传 history = prior', 6.0, 2.65, 2.8, 1.0);
  keyBox(s, '写入最后 tile history', 9.45, 1.9, 2.0, 0.85);
  box(s, '供最终输出使用', 11.75, 1.9, 1.0, 0.85, { fontSize: 11 });
  arrow(s, 2.55, 2.22, 3.0, 2.22); arrow(s, 5.0, 2.1, 6.0, 1.65); arrow(s, 5.0, 2.45, 6.0, 3.05); arrow(s, 8.8, 1.65, 9.45, 2.18); arrow(s, 8.8, 3.05, 9.45, 2.48); arrow(s, 11.45, 2.32, 11.75, 2.32);
  addText(s, '修正后的语义：', 0.75, 4.35, 2.0, 0.3, { fontSize: 15, bold: true });
  bullets(s, ['produced=true：history = prior + new_extrinsic', 'produced=false：history = prior', 'produced=false 只禁止新增译码信息，不禁止保存已有 prior'], 0.95, 4.75, 7.2, 1.25, 15);
  addText(s, '明确保持不变：Unscheduled 仍不执行 EarlyStop/HISO/SISO；produced 仍为 false；tile_out 仍不修改；只修正最后 tile 的 last_tile_history_accum。', 8.4, 4.75, 4.1, 1.25, { fontSize: 13, valign: 'top' });
}

pptx.writeFile({ fileName: 'ofec_editable_sample.pptx' });
