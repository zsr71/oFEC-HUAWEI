const pptxgen = require('pptxgenjs');

const pptx = new pptxgen();
pptx.layout = 'LAYOUT_WIDE';
pptx.author = 'oFEC';
pptx.company = 'oFEC';
pptx.subject = 'Level 5/6 缓冲先进先出共享译码整体方案';
pptx.title = 'oFEC Level 5/6 缓冲先进先出新方案';
pptx.lang = 'zh-CN';
pptx.theme = { headFontFace: 'Microsoft YaHei', bodyFontFace: 'Microsoft YaHei', lang: 'zh-CN' };

const C = { black: '111111', gray: 'F2F2F2', gray2: 'E3E3E3', grid: 'BFC5CC', red: 'C00000', green: '3F7F5F', blue: '536B82', muted: '666666', white: 'FFFFFF', paleRed: 'F9EEEE', paleBlue: 'EEF3F6' };
const FONT = 'Microsoft YaHei';

pptx.defineSlideMaster({
  title: 'MASTER',
  background: { color: C.white },
  objects: [
    { line: { x: 0.75, y: 0.96, w: 11.85, h: 0, line: { color: '222222', width: 1.5 } } },
    { text: { text: 'oFEC · Level 5/6 缓冲先进先出新方案', options: { x: 8.85, y: 7.12, w: 3.7, h: 0.18, fontFace: FONT, fontSize: 8, color: C.muted, align: 'right', margin: 0 } } },
  ],
  slideNumber: { x: 12.62, y: 7.12, color: C.muted, fontFace: FONT, fontSize: 8 },
});

function addText(s, text, x, y, w, h, opts = {}) {
  s.addText(text, { x, y, w, h, fontFace: FONT, fontSize: 14, color: C.black, margin: 0.04, fit: 'shrink', valign: 'mid', breakLine: false, ...opts });
}
function title(s, text, subtitle = '') {
  addText(s, text, 0.75, 0.42, 11.85, 0.4, { fontSize: 22, bold: true });
  if (subtitle) addText(s, subtitle, 0.75, 1.08, 11.7, 0.28, { fontSize: 13.5, color: C.muted });
}
function box(s, text, x, y, w, h, opts = {}) {
  s.addShape(pptx.ShapeType.rect, { x, y, w, h, fill: { color: opts.fill || C.gray }, line: { color: opts.line || '777777', width: opts.lineWidth || 0.8 } });
  addText(s, text, x + 0.08, y + 0.05, w - 0.16, h - 0.1, { align: opts.align || 'center', valign: 'mid', bold: opts.bold || false, fontSize: opts.fontSize || 13, color: opts.color || C.black });
}
function keyBox(s, text, x, y, w, h, opts = {}) { box(s, text, x, y, w, h, { ...opts, line: opts.line || C.red, lineWidth: opts.lineWidth || 1.5, bold: true, fill: opts.fill || C.paleRed }); }
function arrow(s, x1, y1, x2, y2, color = '555555') { s.addShape(pptx.ShapeType.line, { x: x1, y: y1, w: x2 - x1, h: y2 - y1, line: { color, width: 1, endArrowType: 'triangle' } }); }
function bullets(s, items, x, y, w, h, fs = 14) {
  const runs = [];
  items.forEach((item) => { runs.push({ text: item, options: { bullet: { indent: 12 }, hanging: 3, breakLine: true } }); });
  addText(s, runs, x, y, w, h, { fontSize: fs, valign: 'top', paraSpaceAfterPt: 6 });
}
function table(s, rows, x, y, w, h, colW, fs = 12.5, rowH = 0.44) {
  s.addTable(rows, { x, y, w, h, border: { type: 'solid', color: C.grid, pt: 0.6 }, fill: C.white, color: C.black, fontFace: FONT, fontSize: fs, margin: 0.07, valign: 'mid', autoFit: false, rowH, colW, fillHeader: C.gray2, boldHeader: true });
}
function sectionLabel(s, text, x, y, w) { addText(s, text, x, y, w, 0.28, { fontSize: 15, bold: true, color: C.red }); }

// 1. Cover
{
  const s = pptx.addSlide('MASTER');
  s.background = { color: 'EEEEEE' };
  addText(s, 'oFEC 五 / 六级共享译码', 0.85, 1.25, 10.8, 0.42, { fontSize: 20, color: C.muted });
  addText(s, '缓冲先进先出（Buffered FIFO）\n新方案', 0.85, 1.9, 8.8, 1.35, { fontSize: 34, bold: true, breakLine: true, valign: 'top' });
  addText(s, '用跨时刻缓冲吸收低提前停止率造成的瞬时译码需求突发', 0.9, 3.55, 9.0, 0.4, { fontSize: 18, color: C.blue });
  keyBox(s, '64 个码字 / 批次', 0.9, 4.55, 2.2, 0.65, { fontSize: 15 });
  box(s, '先进先出（FIFO）', 3.35, 4.55, 2.2, 0.65, { fontSize: 15 });
  box(s, '有限服务机会', 5.8, 4.55, 2.55, 0.65, { fontSize: 15 });
  addText(s, '2026-08-19  ·  方案评审稿', 0.9, 6.45, 4.5, 0.3, { fontSize: 12.5, color: C.muted });
  // editable visual: burst becomes backlog, then recovers
  const xs = [9.6, 10.35, 11.1, 11.85];
  const heights = [0.6, 1.45, 1.0, 0.45];
  xs.forEach((x, i) => { s.addShape(pptx.ShapeType.rect, { x, y: 5.75 - heights[i], w: 0.52, h: heights[i], fill: { color: i === 1 ? C.red : '777777' }, line: { color: C.white, transparency: 100 } }); });
  addText(s, '需求突发', 9.45, 5.95, 1, 0.25, { fontSize: 11, color: C.muted, align: 'center' });
  addText(s, '恢复', 11.65, 5.95, 1, 0.25, { fontSize: 11, color: C.muted, align: 'center' });
}

// 2. Goal only
{
  const s = pptx.addSlide('MASTER'); title(s, '方案目标', '本页只定义要解决的问题和希望达到的效果');
  keyBox(s, '核心目标', 0.9, 1.7, 1.55, 0.65, { fontSize: 16 });
  addText(s, '把某一时刻突然增加的译码需求，转化为后续时刻可以继续处理的待处理内容。', 2.75, 1.78, 8.8, 0.45, { fontSize: 19, bold: true });
  const goals = [
    ['短时突发可吸收', '需求高的时刻允许暂时积压，不要求当前时刻全部完成。'],
    ['处理顺序不改变', '始终按照批次到达顺序服务，保证先进先出。'],
    ['低负载可追赶', '后续出现低负载时，利用额外完成机会逐步消化积压。'],
    ['边界行为可解释', '缓冲用尽时保持固定输入输出节拍，并明确标记边界顶出。'],
  ];
  goals.forEach((g, i) => { const x = 0.95 + (i % 2) * 6.0; const y = 3.0 + Math.floor(i / 2) * 1.38; box(s, g[0], x, y, 2.0, 0.72, { line: i === 0 ? C.red : C.blue, fill: i === 0 ? C.paleRed : C.paleBlue, bold: true, fontSize: 14 }); addText(s, g[1], x + 2.25, y + 0.05, 3.25, 0.66, { fontSize: 14, valign: 'mid' }); });
  addText(s, '范围说明：本方案讨论整体窗口、批次顺序和跨时刻服务逻辑，不展开具体硬件结构和代码实现。', 0.95, 6.25, 11.0, 0.4, { fontSize: 14, color: C.muted, fill: { color: C.gray }, margin: 0.1 });
}

// 3. Baseline window
{
  const s = pptx.addSlide('MASTER'); title(s, '先看窗口：缓冲区、Level 6 和 Level 5 的相对位置');
  addText(s, '默认示例使用 32 行缓冲区；Level 6 和 Level 5 各占 22 行。这里先只看“逻辑窗口”，不看窗口内部怎么译码。', 0.8, 1.35, 11.3, 0.45, { fontSize: 15, color: C.muted });
  const barX = 1.0, barY = 2.25, barW = 10.8, total = 76;
  const segs = [{ name: '缓冲区\n32 行', n: 32, color: 'D9D9D9' }, { name: 'Level 6\n22 行', n: 22, color: C.blue }, { name: 'Level 5\n22 行', n: 22, color: '8FA6B5' }];
  let cur = barX; segs.forEach((seg) => { const w = barW * seg.n / total; s.addShape(pptx.ShapeType.rect, { x: cur, y: barY, w, h: 0.92, fill: { color: seg.color }, line: { color: C.white, width: 1 } }); addText(s, seg.name, cur + 0.04, barY + 0.1, w - 0.08, 0.65, { fontSize: 14, bold: true, color: seg.color === C.blue ? C.white : C.black, align: 'center' }); cur += w; });
  addText(s, '基准位置：S₀ = 32，当前窗口覆盖最后 44 行', 1.0, 3.55, 6.3, 0.35, { fontSize: 16, bold: true });
  box(s, 'Level 6\n相对固定 22 行', 1.15, 4.25, 2.6, 0.82, { fill: C.paleBlue, line: C.blue, fontSize: 14 });
  box(s, 'Level 5\n相对固定 22 行', 4.15, 4.25, 2.6, 0.82, { fill: 'E8EEF1', line: C.blue, fontSize: 14 });
  arrow(s, 7.1, 4.66, 8.05, 4.66);
  keyBox(s, '窗口起点 S_t\n表示窗口向缓冲区回退了多少', 8.05, 4.0, 3.2, 1.32, { fontSize: 15 });
  addText(s, '窗口移动时，Level 5 和 Level 6 一起移动；两者之间的相对布局不改变。', 1.0, 5.85, 11.0, 0.42, { fontSize: 15, color: C.blue, fill: { color: C.paleBlue }, margin: 0.1 });
}

// 4. Physical movement per time step
{
  const s = pptx.addSlide('MASTER'); title(s, '每个时刻的物理移动：上面出去两行，下面进来两行');
  addText(s, '无论当前译码是否完成批次，物理阵列每个时刻都执行同一个动作。', 0.8, 1.35, 10.8, 0.35, { fontSize: 15, color: C.muted });
  const rows = [
    ['时刻开始', '顶部 [0,2) 的两行准备输出', '底部 [74,76) 保留当前最新输入'],
    ['固定移动', '顶部两行输出；其余内容整体上移两行', '原来的 [2,76) → [0,74)'],
    ['新输入到达', '底部空出的两行写入新批次', '新批次两行进入 [74,76)'],
  ];
  rows.forEach((r, i) => { const y = 2.05 + i * 1.05; keyBox(s, r[0], 0.9, y, 1.55, 0.62, { fontSize: 14, fill: i === 1 ? C.paleRed : C.white }); box(s, r[1], 2.8, y, 3.55, 0.62, { fontSize: 13 }); arrow(s, 6.55, y + 0.31, 7.15, y + 0.31); box(s, r[2], 7.25, y, 4.45, 0.62, { fill: C.paleBlue, line: C.blue, fontSize: 13 }); });
  addText(s, '以连续时刻为例：', 0.9, 5.35, 2.0, 0.3, { fontSize: 16, bold: true });
  box(s, 't=0：B0 两行进入底部', 1.05, 5.85, 2.6, 0.62, { fontSize: 13 }); arrow(s, 3.8, 6.16, 4.35, 6.16); box(s, 't=1：B1 进入底部，B0 上移两行', 4.4, 5.85, 3.3, 0.62, { fontSize: 13 }); arrow(s, 7.85, 6.16, 8.4, 6.16); box(s, 't=2：B2 进入底部，B0 再上移两行', 8.45, 5.85, 3.3, 0.62, { fontSize: 13 });
}

// 5. Window delay and recovery example
{
  const s = pptx.addSlide('MASTER'); title(s, '窗口如何延后，又如何提前恢复？', '直接使用方案文档中的 t=0 到 t=4 逻辑算例');
  table(s, [
    [{ text: '时刻', options: { bold: true } }, { text: '队首批次', options: { bold: true } }, { text: '窗口起点', options: { bold: true } }, { text: '窗口变化原因', options: { bold: true } }, { text: '窗口状态', options: { bold: true } }],
    ['t=0', 'B0', '32 → 30', 'B0 需求突增，本时刻没有完成批次', '向缓冲区回退 2 行'],
    ['t=1', 'B0', '30 → 28', 'B0 仍有待处理内容，继续优先服务 B0', '再次回退 2 行'],
    ['t=2', 'B0 完成', '28 → 28', '完成 1 个批次，刚好抵消固定推进', '窗口保持不动'],
    ['t=3', 'B1、B2 快速完成', '28 → 30', '低负载，连续完成 2 个批次', '向基准位置恢复 2 行'],
    ['t=4', 'B3、B4 完成', '30 → 32', '继续出现可快速完成的批次', '恢复到基准位置'],
  ], 0.7, 1.55, 12.0, 3.45, [0.8, 2.15, 1.35, 4.35, 3.25], 11.5, 0.52);
  addText(s, '核心规律', 0.85, 5.35, 1.5, 0.3, { fontSize: 16, bold: true, color: C.red });
  bullets(s, ['需求高、批次完成少：窗口向缓冲区回退，给待处理内容留下空间。', '一个时刻完成 1 个批次：只抵消固定的两行物理推进，窗口不恢复。', '一个时刻额外完成多个批次：窗口才向基准位置提前恢复。'], 1.0, 5.75, 10.9, 0.9, 14);
}

// 6. Window transition intuition
{
  const s = pptx.addSlide('MASTER'); title(s, '窗口移动的直观规则');
  keyBox(s, '每个时刻先固定上移两行，再看内部完成了多少批次', 2.35, 1.55, 8.15, 0.78, { fontSize: 17 });
  const cases = [
    ['完成 0 个批次', '窗口向后退 2 行', '留下空间，继续跟随队首待处理内容', C.red],
    ['完成 1 个批次', '窗口位置不变', '批次完成带来的推进，刚好抵消固定上移', C.blue],
    ['完成 2 个及以上', '窗口向前恢复', '多出来的完成量用于消化之前的积压', C.green],
  ];
  cases.forEach((c, i) => { const x = 0.9 + i * 4.05; box(s, c[0], x, 3.0, 1.55, 0.65, { line: c[3], fill: C.white, bold: true, fontSize: 13 }); box(s, c[1], x + 1.7, 3.0, 1.8, 0.65, { fill: i === 0 ? C.paleRed : i === 1 ? C.paleBlue : 'EEF5EF', line: c[3], fontSize: 13 }); addText(s, c[2], x, 3.9, 3.55, 0.7, { fontSize: 13, valign: 'top', color: C.muted }); });
  addText(s, '这一步只描述“窗口怎么动”，还没有进入单个批次内部的译码细节。', 0.9, 5.35, 11.4, 0.45, { fontSize: 15, color: C.blue, fill: { color: C.paleBlue }, margin: 0.1 });
  addText(s, '窗口恢复不等于对外一次输出多个批次：物理输入输出节拍始终保持每个时刻两行。', 0.9, 6.02, 11.4, 0.4, { fontSize: 14, bold: true });
}

// 7. Service flow
{
  const s = pptx.addSlide('MASTER'); title(s, '然后再看：一个时刻内部如何服务队首批次');
  const stages = [
    ['队首批次', '先进先出'], ['跳过已完成内容', '不重复处理'], ['原共享译码流程', '提前停止 / 分类 / 调度'], ['结果分流', '完成 / 待处理'], ['窗口更新', '完成数决定方向'],
  ];
  stages.forEach((st, i) => { const x = 0.85 + i * 2.42; box(s, st[0] + '\n' + st[1], x, 1.75, 1.78, 0.98, { line: i === 2 ? C.red : '777777', fill: i === 2 ? C.paleRed : C.gray, bold: i === 2, fontSize: 13 }); if (i < stages.length - 1) arrow(s, x + 1.82, 2.24, x + 2.35, 2.24); });
  sectionLabel(s, '普通批次路径', 0.85, 3.35, 2.2);
  bullets(s, ['本时刻只使用当前批次的有限服务机会。', '已经完成的码字跳过；未完成的码字按原共享译码流程继续处理。', '本时刻仍未完成的内容保留到下一时刻，不丢失、不重复。'], 1.0, 3.75, 5.5, 1.2, 14);
  sectionLabel(s, '全部提前停止的快速路径', 7.0, 3.35, 3.4);
  keyBox(s, '64 个码字全部提前停止', 7.1, 3.78, 2.55, 0.72, { fontSize: 14 });
  arrow(s, 9.8, 4.14, 10.35, 4.14);
  box(s, '立即完成\n计入完成数', 10.35, 3.78, 1.55, 0.72, { fontSize: 14 });
  addText(s, '快速路径不占普通共享译码机会；同一时刻可以继续检查下一个队首批次。', 7.1, 4.85, 5.0, 0.58, { fontSize: 14, color: C.muted, valign: 'top' });
}

// 6. Cross-time logic
{
  const s = pptx.addSlide('MASTER'); title(s, '跨时刻处理逻辑：完成的内容退出，未完成的内容留下');
  table(s, [
    [{ text: '上一时刻的结果', options: { bold: true } }, { text: '下一时刻怎么处理', options: { bold: true } }, { text: '方案含义', options: { bold: true } }],
    ['已经完成的码字', '直接跳过', '不再重复译码，也不再重复写回'],
    ['尚未完成的码字', '继续参与本批次服务', '保留原来的处理状态，等待后续机会'],
    ['整批次全部提前停止', '立即完成并转向下一批次', '利用低负载时刻快速消化队列'],
  ], 0.85, 1.55, 11.4, 2.35, [2.55, 3.55, 5.3], 13, 0.55);
  addText(s, '原五六级共享译码流程不变', 0.85, 4.35, 3.5, 0.3, { fontSize: 16, bold: true, color: C.red });
  const rules = [['提前停止判断', '先判断哪些码字可以直接结束'], ['难度分类', '对剩余码字判断适合的译码路径'], ['有限服务机会', '按原规则服务当前队首批次']];
  rules.forEach((r, i) => { box(s, r[0], 1.0 + i * 3.82, 4.9, 1.75, 0.65, { fontSize: 13, bold: true }); arrow(s, 2.8 + i * 3.82, 5.22, 3.15 + i * 3.82, 5.22); box(s, r[1], 3.15 + i * 3.82, 4.9, 2.0, 0.65, { fill: C.paleBlue, line: C.blue, fontSize: 13 }); });
  addText(s, '新增的只是“跨时刻保留 + 已完成跳过”这一层逻辑，不改变原有译码算法本身。', 0.85, 6.15, 11.3, 0.4, { fontSize: 14, color: C.blue, fill: { color: C.paleBlue }, margin: 0.1 });
}

// 8. Formula
{
  const s = pptx.addSlide('MASTER'); title(s, '窗口状态转移：完成多少，决定回退还是恢复');
  keyBox(s, 'S_{t+1} = clip(S_t - 2 + 2·C_t, 0, R_buf)', 2.25, 1.55, 8.1, 0.85, { fontSize: 21, fill: C.paleRed });
  const cases = [
    ['完成数 = 0', '窗口回退', '没有批次完成，继续跟随待处理队首', C.red],
    ['完成数 = 1', '窗口保持', '一个批次完成，抵消固定推进', C.blue],
    ['完成数 ≥ 2', '窗口恢复', '额外完成批次，窗口向基准位置移动', C.green],
  ];
  cases.forEach((c, i) => { const x = 0.9 + i * 4.05; box(s, c[0], x, 3.0, 1.25, 0.62, { line: c[3], fill: C.white, bold: true, fontSize: 14 }); box(s, c[1], x + 1.4, 3.0, 2.15, 0.62, { fill: i === 0 ? C.paleRed : i === 1 ? C.paleBlue : 'EEF5EF', line: c[3], fontSize: 13 }); addText(s, c[2], x, 3.85, 3.55, 0.7, { fontSize: 13, valign: 'top', color: C.muted }); });
  addText(s, '完成数包含普通完成和“全部提前停止”的快速完成，不等于对外接口一次输出多少批次。', 0.9, 5.25, 11.4, 0.42, { fontSize: 14, color: C.blue, fill: { color: C.paleBlue }, margin: 0.1 });
  addText(s, '对外输入输出节拍保持固定；方案只改变内部如何消化积压。', 0.9, 5.95, 11.4, 0.42, { fontSize: 14, bold: true });
}

// 9. Boundary
{
  const s = pptx.addSlide('MASTER'); title(s, '缓冲用尽时怎么办？', '这是有限缓冲方案必须明确的边界，不是正常的完成路径');
  keyBox(s, '窗口已经退到最前端，队首批次仍有待处理内容', 2.45, 1.55, 6.9, 0.78, { fontSize: 17 });
  const steps = [['1', '输入输出节拍照常进行，不暂停系统'], ['2', '队首仍未完成的内容被顶出，标记为边界顶出（forced-evicted）'], ['3', '顶出内容保留当前结果，不再等待新的译码机会'], ['4', '边界顶出不算正常完成，但结果仍纳入统一质量统计']];
  steps.forEach((st, i) => { keyBox(s, st[0], 1.0, 2.8 + i * 0.56, 0.36, 0.34, { fontSize: 11, fill: C.white }); addText(s, st[1], 1.55, 2.8 + i * 0.56, 9.8, 0.34, { fontSize: 14, color: i === 2 ? C.red : C.black, bold: i === 2 }); });
  box(s, '边界顶出\n是保持固定节拍的最后保护\n不是正常完成，也不是提前停止完成', 8.65, 4.9, 3.0, 1.15, { fill: C.paleRed, line: C.red, fontSize: 13, bold: true });
  addText(s, '缓冲深度越大，能够吸收的短时缺口越大；但缓冲不是无限的。', 1.0, 5.95, 7.2, 0.35, { fontSize: 14, color: C.muted });
}

// 10. Boundary timeline example
{
  const s = pptx.addSlide('MASTER'); title(s, '缓冲用尽的连续示例：B0、B1、B2 如何依次获得服务', '采用较小的示例缓冲区 R_buf=4 行；表格与方案文档的逐时刻案例采用相同呈现方式');
  table(s, [
    [{ text: '时刻', options: { bold: true } }, { text: '新到达', options: { bold: true } }, { text: 'FIFO 队列', options: { bold: true } }, { text: '当前队首', options: { bold: true } }, { text: 'B0 状态', options: { bold: true } }, { text: 'B1 状态', options: { bold: true } }, { text: 'B2 状态', options: { bold: true } }, { text: '窗口 / 本时刻结果', options: { bold: true } }],
    ['t=0', 'B0', '[B0]', 'B0', '首次服务后仍未完成', '—', '—', 'S=4 → 2；窗口回退'],
    ['t=1', 'B1', '[B0, B1]', 'B0', '再次服务后仍未完成', '等待', '—', 'S=2 → 0；窗口到边界'],
    ['t=2', 'B2', '[B0, B1, B2]', 'B0', '第三次服务后仍未完成\n→ 边界顶出', '等待', '等待', 'S=0；B0 顶出，不计正常完成'],
    ['t=3', 'B3', '[B1, B2, B3]', 'B1', '已顶出', '获得自己的首次服务\n（示例中完成）', '等待', 'S=0；B1 不会被 B0 剥夺机会'],
    ['t=4', 'B4', '[B2, B3, B4]', 'B2', '已顶出', '已完成', '获得自己的首次服务\n（示例中完成）', 'S=0；FIFO 继续向后服务'],
  ], 0.42, 1.42, 12.45, 3.85, [0.6, 0.82, 1.65, 0.9, 1.85, 1.8, 1.9, 2.93], 10.6, 0.55);
  addText(s, '读表重点', 0.75, 5.58, 1.35, 0.3, { fontSize: 16, bold: true, color: C.red });
  bullets(s, ['B0 在 t=0、t=1、t=2 始终是队首，因此 B1、B2 只能等待，不能越过 B0。', 't=2 触及 S=0 后，B0 仍先获得一次普通服务；只有该次服务后仍未完成，才执行边界顶出。', 'B0 顶出后，B1 在 t=3、B2 在 t=4 依次获得各自的首次服务机会；等待是延迟，不是取消服务机会。'], 0.9, 5.95, 11.5, 0.8, 13.4);
}

// 11. Invariants and non-goals
{
  const s = pptx.addSlide('MASTER'); title(s, '方案边界：必须保持的整体逻辑');
  table(s, [
    [{ text: '不变量', options: { bold: true } }, { text: '含义', options: { bold: true } }],
    ['先进先出（FIFO）', 'B0 → B1 → B2；队首未完成时，后续批次不越过'],
    ['一次完成一次处理', '一个码字完成后退出，不重复处理、不重复写回'],
    ['完成内容退出', '后续时刻只处理仍未完成的内容'],
    ['窗口相对关系稳定', 'Level 5 / Level 6 的逻辑窗口大小和相对位置不变'],
    ['输入输出节拍稳定', '每个时刻仍按固定节拍接收和输出'],
    ['完成统计清晰', '正常完成和快速完成计入完成数；边界顶出单独统计'],
  ], 0.85, 1.45, 11.5, 3.75, [2.0, 9.5], 13, 0.55);
  addText(s, '明确不采用', 0.85, 5.55, 1.8, 0.3, { fontSize: 16, bold: true, color: C.red });
  bullets(s, ['后续批次越过未完成的先进先出队首。', '把一个批次剩余的服务机会直接转给另一个普通批次。', '根据未来批次的提前停止率提前做预测性调度。', '在等待期间改变原有批次的处理顺序。'], 1.0, 5.95, 10.8, 0.8, 13.5);
}

// 12. Verification and close
{
  const s = pptx.addSlide('MASTER'); title(s, '如何验证方案是否有效？');
  const cols = [
    ['功能正确性', ['每个码字最终只有一种结果：正常完成、提前停止完成或边界顶出', '已完成内容不重复处理', '先进先出顺序和待处理延后符合定义']],
    ['窗口与吞吐', ['高负载时窗口会回退，低负载时能够恢复', '窗口不会越过缓冲边界', '输入输出节拍始终保持稳定']],
    ['压力与质量', ['观察积压峰值、最老等待时间和恢复时间', '单独统计边界顶出数量', '统一比较整体译码质量']],
  ];
  cols.forEach((c, i) => { const x = 0.85 + i * 4.0; keyBox(s, c[0], x, 1.55, 3.25, 0.62, { fontSize: 15, fill: i === 1 ? C.paleBlue : C.paleRed, line: i === 1 ? C.blue : C.red }); bullets(s, c[1], x + 0.15, 2.45, 3.0, 1.5, 13.2); });
  addText(s, '方案结论', 0.85, 4.55, 1.5, 0.3, { fontSize: 16, bold: true, color: C.red });
  addText(s, '缓冲先进先出方案把“某一时刻处理不完的需求”变成“后续时刻继续处理的待处理内容”，并利用低负载时刻的快速完成路径逐步消化积压。它解决的是短时波动，不是长期平均能力不足。', 0.95, 5.0, 11.0, 0.85, { fontSize: 17, bold: true, color: C.black, fill: { color: C.gray }, margin: 0.14, valign: 'mid' });
  addText(s, '建议评审重点：整体流程是否闭环、先进先出是否清晰、积压能否恢复、边界是否可控。', 0.95, 6.18, 11.0, 0.38, { fontSize: 14, color: C.blue });
}

pptx.writeFile({ fileName: __dirname + '/ofec_level56_buffered_fifo_new_scheme.pptx' });
