import {DualWindowSimulator, SCENARIOS, runSelfTests} from './simulator.js';

const byId = (id) => document.getElementById(id);
const elements = {
  scenario: byId('scenarioSelect'),
  rbuf: byId('rbufInput'),
  reset: byId('resetButton'),
  previous: byId('previousButton'),
  next: byId('nextButton'),
  nextTime: byId('nextTimeButton'),
  play: byId('playButton'),
  speed: byId('speedInput'),
  speedValue: byId('speedValue'),
  memory: byId('memoryGrid'),
  memorySize: byId('memorySize'),
  windowL5: byId('windowL5'),
  windowL6: byId('windowL6'),
  timeTimeline: byId('timeTimeline'),
  eventTimeline: byId('eventTimeline'),
  scenarioDescription: byId('scenarioDescription'),
  entryBody: byId('entryTableBody'),
  frameCounter: byId('frameCounter'),
  budgetCounter: byId('budgetCounter'),
  roundSummary: byId('roundSummary'),
  eventTitle: byId('eventTitle'),
  eventDescription: byId('eventDescription'),
  eventPhase: byId('eventPhase'),
  positionState: byId('positionState'),
  headState: byId('headState'),
  recoveryState: byId('recoveryState'),
  forcedState: byId('forcedState'),
  assertionList: byId('assertionList'),
  assertionSummary: byId('assertionSummary'),
  selfTestBadge: byId('selfTestBadge'),
  metricTime: byId('metricTime'),
  metricEvent: byId('metricEvent'),
  metricL5: byId('metricL5'),
  metricL6: byId('metricL6'),
  metricEntries: byId('metricEntries'),
  metricCompletion: byId('metricCompletion'),
  positionStartLabel: byId('positionStartLabel'),
  positionStartL5: byId('positionStartL5'),
  positionStartL6: byId('positionStartL6'),
  positionRoundL5: byId('positionRoundL5'),
  positionRoundL6: byId('positionRoundL6'),
  positionCommitBox: byId('positionCommitBox'),
  positionCommitLabel: byId('positionCommitLabel'),
  positionCommitL5: byId('positionCommitL5'),
  positionCommitL6: byId('positionCommitL6'),
};

let simulation;
let trace;
let frameIndex = 0;
let timer = null;

function batchText(id) {
  if (id === null || id === undefined) return '—';
  return id >= 0 ? `B${id}` : `B−${Math.abs(id)}`;
}

function rangeText([start, end]) {
  return `[${start},${end})`;
}

function escapeHtml(value) {
  return String(value)
    .replaceAll('&', '&amp;')
    .replaceAll('<', '&lt;')
    .replaceAll('>', '&gt;')
    .replaceAll('"', '&quot;')
    .replaceAll("'", '&#039;');
}

function currentSnapshot() {
  return trace.snapshots[frameIndex];
}

function stopPlaying() {
  if (timer !== null) window.clearInterval(timer);
  timer = null;
  elements.play.textContent = '播放';
}

function startPlaying() {
  stopPlaying();
  elements.play.textContent = '暂停';
  timer = window.setInterval(() => {
    if (frameIndex >= trace.snapshots.length - 1) {
      stopPlaying();
      return;
    }
    frameIndex += 1;
    render();
  }, Number(elements.speed.value));
}

function buildSimulation({preserveScenario = true} = {}) {
  stopPlaying();
  const scenarioId = preserveScenario ? elements.scenario.value : SCENARIOS[0].id;
  const scenario = SCENARIOS.find((item) => item.id === scenarioId) || SCENARIOS[0];
  const rbuf = Math.max(0, Math.floor(Number(elements.rbuf.value || scenario.defaultRbuf) / 2) * 2);
  elements.rbuf.value = String(rbuf);
  simulation = new DualWindowSimulator({scenarioId: scenario.id, rbuf});
  trace = simulation.simulate(scenario.times);
  elements.scenarioDescription.textContent = scenario.description;
  frameIndex = 0;
  render();
}

function regionForRow(snapshot, row) {
  const inD5 = row >= snapshot.d5[0] && row < snapshot.d5[1];
  const inD6 = row >= snapshot.d6[0] && row < snapshot.d6[1];
  if (inD5) return 'd5';
  if (inD6) return 'd6';
  if (row >= snapshot.s5 && row < snapshot.s5 + 22) return 'l5';
  if (row >= snapshot.s6 && row < snapshot.s6 + 22) return 'l6';
  return 'buffer';
}

function renderMemory(snapshot) {
  const rowHeight = Math.max(7, Math.min(14, Math.floor(670 / (snapshot.rmem + 1))));
  elements.memory.style.setProperty('--row-height', `${rowHeight}px`);
  elements.memorySize.textContent = `${snapshot.rmem} rows × 8 blocks`;
  elements.windowL5.textContent = `当前轮次 W5 [${snapshot.s5},${snapshot.s5 + 22})`;
  elements.windowL6.textContent = `当前轮次 W6 [${snapshot.s6},${snapshot.s6 + 22})`;

  const heads = ['row', 'B', 'half', ...Array.from({length: 8}, (_, index) => `b${index}`)];
  let html = heads.map((head) => `<div class="grid-head" role="columnheader">${head}</div>`).join('');
  for (let row = 0; row < snapshot.rmem; row += 1) {
    const data = snapshot.memory[row];
    const region = regionForRow(snapshot, row);
    const input = data.batch === snapshot.incomingBatch && row >= snapshot.rmem - 2;
    const output = row < 2;
    const extra = input ? ' input' : (output ? ' output' : '');
    const title = `row ${row} / ${batchText(data.batch)} / half ${data.half}`;
    html += `<div class="row-label ${region}${extra}" title="${title}">${row}</div>`;
    html += `<div class="batch-label ${region}${extra}" title="${title}">${batchText(data.batch)}</div>`;
    html += `<div class="half-label ${region}${extra}" title="${title}">r${data.half}</div>`;
    for (let block = 0; block < 8; block += 1) {
      const classes = ['memory-cell', region === 'buffer' ? (data.batch % 2 === 0 ? 'base-a' : 'base-b') : region];
      if (row === snapshot.s5 || row === snapshot.s6) classes.push('window-top');
      if (row === snapshot.s5 + 21 || row === snapshot.s6 + 21) classes.push('window-bottom');
      if (input) classes.push('input');
      else if (output) classes.push('output');
      html += `<div class="${classes.join(' ')}" role="gridcell" title="${title}, block ${block}"></div>`;
    }
  }
  elements.memory.innerHTML = html;
}

function renderMetrics(snapshot) {
  elements.metricTime.textContent = `t=${snapshot.t}`;
  elements.metricEvent.textContent = snapshot.event.title;
  elements.metricL5.textContent = `${batchText(snapshot.h5)} · [${snapshot.s5},${snapshot.s5 + 22})`;
  elements.metricL6.textContent = `${batchText(snapshot.h6)} · [${snapshot.s6},${snapshot.s6 + 22})`;
  elements.metricEntries.textContent = `${snapshot.entries.length} / 8`;
  elements.metricCompletion.textContent = `C5=${snapshot.c5} · C6=${snapshot.c6}`;
}

function renderPositionStrip(snapshot) {
  const start5 = snapshot.timeStartS5;
  const start6 = snapshot.timeStartS6;
  elements.positionStartLabel.textContent = `t=${snapshot.t} 开始位置`;
  elements.positionStartL6.textContent = `L6 S6_t=${start6} · W6 [${start6},${start6 + 22})`;
  elements.positionStartL5.textContent = `L5 S5_t=${start5} · W5 [${start5},${start5 + 22})`;
  elements.positionRoundL6.textContent = `L6 S6^(r)=${snapshot.s6} · W6 [${snapshot.s6},${snapshot.s6 + 22})`;
  elements.positionRoundL5.textContent = `L5 S5^(r)=${snapshot.s5} · W5 [${snapshot.s5},${snapshot.s5 + 22})`;

  const resolved = snapshot.nextS5 !== null && snapshot.nextS6 !== null;
  elements.positionCommitBox.classList.toggle('resolved', resolved);
  elements.positionCommitLabel.textContent = resolved ? `t=${snapshot.t + 1} 提交位置` : `t=${snapshot.t + 1} 提交位置`;
  elements.positionCommitL6.textContent = resolved
    ? `L6 S6_${snapshot.t + 1}=${snapshot.nextS6} · W6 [${snapshot.nextS6},${snapshot.nextS6 + 22})`
    : 'L6 待时刻边界联合裁决';
  elements.positionCommitL5.textContent = resolved
    ? `L5 S5_${snapshot.t + 1}=${snapshot.nextS5} · W5 [${snapshot.nextS5},${snapshot.nextS5 + 22})`
    : 'L5 待时刻边界联合裁决';
}

function renderTimeTimeline(snapshot) {
  const available = [...trace.summaries];
  if (!available.some((item) => item.t === snapshot.t)) {
    available.push({
      t: snapshot.t,
      inputBatch: snapshot.incomingBatch,
      h5Before: snapshot.h5,
      h6Before: snapshot.h6,
      s5Before: snapshot.committedS5,
      s6Before: snapshot.committedS6,
      s5After: snapshot.s5,
      s6After: snapshot.s6,
      c5: snapshot.c5,
      c6: snapshot.c6,
      entries: snapshot.entries.length,
      rounds: snapshot.rounds,
      finalSnapshotIndex: frameIndex,
    });
  }
  elements.timeTimeline.innerHTML = available.map((item) => `
    <button class="time-node${item.t === snapshot.t ? ' active' : ''}" type="button" data-time="${item.t}" title="查看 t=${item.t}">
      <strong>t=${item.t} · 输入 ${batchText(item.inputBatch)}</strong>
      <span>L5 ${item.s5Before}→${item.s5After} · L6 ${item.s6Before}→${item.s6After}</span>
      <span class="completion">C5/C6=${item.c5}/${item.c6} · entry=${item.entries} · ${item.rounds}轮</span>
    </button>`).join('');
  elements.timeTimeline.querySelectorAll('.time-node').forEach((button) => {
    button.addEventListener('click', () => {
      const time = Number(button.dataset.time);
      const summary = trace.summaries.find((item) => item.t === time);
      if (summary) frameIndex = summary.finalSnapshotIndex;
      else frameIndex = trace.snapshots.findIndex((item) => item.t === time);
      render();
    });
  });
  elements.timeTimeline.querySelector('.time-node.active')?.scrollIntoView({block: 'nearest', inline: 'center'});
}

function renderEventTimeline(snapshot) {
  const indices = trace.snapshots
    .map((item, index) => ({item, index}))
    .filter(({item}) => item.t === snapshot.t);
  elements.eventTimeline.innerHTML = indices.map(({item, index}) => `
    <button class="event-node${index === frameIndex ? ' active' : ''}" type="button" data-index="${index}">
      ${escapeHtml(item.event.title)}
    </button>`).join('');
  elements.eventTimeline.querySelectorAll('.event-node').forEach((button) => {
    button.addEventListener('click', () => {
      frameIndex = Number(button.dataset.index);
      render();
    });
  });
  elements.eventTimeline.querySelector('.event-node.active')?.scrollIntoView({block: 'nearest', inline: 'center'});
}

function renderEntries(snapshot) {
  const rows = [];
  for (let slot = 0; slot < 8; slot += 1) {
    const entry = snapshot.entries.find((item) => item.slot === slot);
    if (!entry) {
      rows.push(`<tr class="empty"><td>E${slot}</td><td>—</td><td>未使用</td><td>—</td><td>—</td><td>—</td></tr>`);
      continue;
    }
    rows.push(`<tr class="l${entry.level}"><td>E${entry.slot}</td><td>R${entry.round}</td><td>L${entry.level} / ${batchText(entry.batch)}</td><td>G${entry.group}</td><td>code ${entry.hiso}</td><td>code ${entry.siso}</td></tr>`);
  }
  elements.entryBody.innerHTML = rows.join('');
  elements.budgetCounter.textContent = `剩余 ${snapshot.remainingEntries}`;
  elements.roundSummary.textContent = snapshot.rounds > 0
    ? `已执行 ${snapshot.rounds} 轮；各轮共享同一份 8-entry 总预算`
    : '尚未执行普通轮次';
}

function renderEventDetail(snapshot) {
  elements.eventTitle.textContent = snapshot.event.title;
  elements.eventDescription.textContent = snapshot.event.detail;
  elements.eventPhase.textContent = snapshot.event.phase.toUpperCase();
  elements.positionState.textContent = `S6=${snapshot.s6}，S5=${snapshot.s5}`;
  elements.headState.textContent = `L5(${batchText(snapshot.h5)})，L6(${batchText(snapshot.h6)})`;
  elements.recoveryState.textContent = `L5=${snapshot.recovery[5]} 行，L6=${snapshot.recovery[6]} 行`;
  const forced = [
    ...snapshot.forced5.map((batch) => `L5(${batchText(batch)})`),
    ...snapshot.forced6.map((batch) => `L6(${batchText(batch)})`),
  ];
  elements.forcedState.textContent = forced.length > 0 ? forced.join('、') : '无';
}

function renderAssertions(snapshot) {
  const failures = snapshot.assertions.filter((item) => !item.pass);
  elements.assertionSummary.textContent = failures.length === 0 ? '当前帧全部通过' : `${failures.length} 项失败`;
  elements.assertionList.innerHTML = snapshot.assertions.map((item) => `
    <div class="assertion-item${item.pass ? '' : ' fail'}">
      <div class="mark">${item.pass ? '✓' : '×'}</div>
      <div><b>${escapeHtml(item.name)}</b><span>${escapeHtml(item.detail)}</span></div>
    </div>`).join('');
}

function render() {
  const snapshot = currentSnapshot();
  renderMetrics(snapshot);
  renderPositionStrip(snapshot);
  renderMemory(snapshot);
  renderTimeTimeline(snapshot);
  renderEventTimeline(snapshot);
  renderEntries(snapshot);
  renderEventDetail(snapshot);
  renderAssertions(snapshot);
  elements.frameCounter.textContent = `${frameIndex + 1} / ${trace.snapshots.length}`;
  elements.previous.disabled = frameIndex === 0;
  elements.next.disabled = frameIndex === trace.snapshots.length - 1;
}

function jumpToNextTime() {
  const currentT = currentSnapshot().t;
  const nextIndex = trace.snapshots.findIndex((snapshot, index) => index > frameIndex && snapshot.t > currentT);
  frameIndex = nextIndex >= 0 ? nextIndex : trace.snapshots.length - 1;
  render();
}

function initializeOptions() {
  elements.scenario.innerHTML = SCENARIOS
    .map((scenario) => `<option value="${scenario.id}">${escapeHtml(scenario.name)}</option>`)
    .join('');
  elements.scenario.value = SCENARIOS[0].id;
  elements.rbuf.value = String(SCENARIOS[0].defaultRbuf);
}

function renderSelfTests() {
  const tests = runSelfTests();
  const failed = tests.filter((test) => !test.pass);
  elements.selfTestBadge.textContent = failed.length === 0 ? `${tests.length} 项自检通过` : `${failed.length} 项自检失败`;
  elements.selfTestBadge.classList.add(failed.length === 0 ? 'pass' : 'fail');
  if (failed.length > 0) {
    elements.selfTestBadge.title = failed.map((test) => `${test.name}: ${test.detail}`).join('\n');
  }
}

elements.scenario.addEventListener('change', () => {
  const scenario = SCENARIOS.find((item) => item.id === elements.scenario.value);
  elements.rbuf.value = String(scenario.defaultRbuf);
  buildSimulation();
});
elements.rbuf.addEventListener('change', () => buildSimulation());
elements.reset.addEventListener('click', () => buildSimulation());
elements.previous.addEventListener('click', () => {
  stopPlaying();
  frameIndex = Math.max(0, frameIndex - 1);
  render();
});
elements.next.addEventListener('click', () => {
  stopPlaying();
  frameIndex = Math.min(trace.snapshots.length - 1, frameIndex + 1);
  render();
});
elements.nextTime.addEventListener('click', () => {
  stopPlaying();
  jumpToNextTime();
});
elements.play.addEventListener('click', () => {
  if (timer !== null) stopPlaying();
  else {
    if (frameIndex === trace.snapshots.length - 1) frameIndex = 0;
    startPlaying();
  }
});
elements.speed.addEventListener('input', () => {
  elements.speedValue.textContent = `${(Number(elements.speed.value) / 1000).toFixed(2)} s`;
  if (timer !== null) startPlaying();
});

initializeOptions();
buildSimulation();
renderSelfTests();
