const WINDOW_ROWS = 22;
const OBJECT_ROWS = 2;
const ENTRY_BUDGET = 8;
const BLOCK_COLUMNS = 8;

const clone = (value) => JSON.parse(JSON.stringify(value));

function evenRows(value) {
  const number = Number.isFinite(Number(value)) ? Number(value) : 32;
  return Math.max(0, Math.min(128, Math.floor(number / 2) * 2));
}

function baselinePositions(rbuf) {
  return {s6: rbuf, s5: rbuf + WINDOW_ROWS};
}

function floatedPositions(rbuf) {
  const slackBelow = Math.min(6, rbuf);
  return {
    s6: rbuf - slackBelow,
    s5: rbuf + WINDOW_ROWS - slackBelow,
  };
}

function topPositions() {
  return {s6: 0, s5: WINDOW_ROWS};
}

function ordinary(demand, note = '') {
  return {kind: 'ordinary', demand, note};
}

function fullEarlyStop(note = '') {
  return {kind: 'full-es', demand: 0, note};
}

export const SCENARIOS = [
  {
    id: 'steady',
    name: '稳态：两级同步完成',
    description: '两级各需要 4 个 group-entry。一个 t 内正好用满 8 个 entry，两个窗口在物理推进后保持基准位置。',
    defaultRbuf: 32,
    times: 5,
    initialPositions: baselinePositions,
    profile: ({level}) => ordinary(4, `Level ${level} 稳态负载`),
  },
  {
    id: 'l6-burst',
    name: 'Level 6 burst：Level 5 独立前进',
    description: '首个 Level 6 对象需要 10 个 group-entry，跨 t 保留 pending；Level 5 负载较轻，可以独立完成和移动。',
    defaultRbuf: 32,
    times: 6,
    initialPositions: baselinePositions,
    profile: ({level, offset}) => {
      if (level === 6 && offset === 0) return ordinary(10, 'Level 6 burst');
      return ordinary(level === 5 ? 2 : 4, '普通负载');
    },
  },
  {
    id: 'both-pending',
    name: '相邻且双 pending：同步向顶部移动',
    description: '两个窗口初始相邻，两级需求都超过本 t 预算。边界处两个窗口同步向顶部移动 2 行，不触发 Level 5 强制弹出。',
    defaultRbuf: 32,
    times: 5,
    initialPositions: baselinePositions,
    profile: () => ordinary(16, '持续高负载'),
  },
  {
    id: 'deferred-recovery',
    name: '恢复受阻：保留窗口推进欠账',
    description: 'Level 6 全 EarlyStop 完成，但因窗口相邻暂时不能换入下一对象；推进欠账保留，待 Level 5 完成并腾出空间后再恢复。',
    defaultRbuf: 32,
    times: 5,
    initialPositions: () => ({s6: 2, s5: 24}),
    profile: ({level, offset}) => {
      if (level === 6 && offset === 0) return fullEarlyStop('恢复会暂时受 Level 5 阻挡');
      if (level === 5 && offset === 0) return ordinary(16, '两次 t 才能完成');
      return ordinary(4, '后续普通对象');
    },
  },
  {
    id: 'top-boundary',
    name: 'Level 6 顶部边界 forced-evicted',
    description: 'Level 6 位于 SRAM 顶部并在本 t 服务后仍有 pending。只终止 Level 6 阶段等待，物理数据继续留在 SRAM 中流动。',
    defaultRbuf: 32,
    times: 4,
    initialPositions: topPositions,
    profile: ({level, offset}) => {
      if (level === 5 && offset === 0) return fullEarlyStop('为 Level 5 腾出一组位置');
      if (level === 6 && offset === 0) return ordinary(16, '顶部持续 pending');
      return ordinary(4, '后续普通对象');
    },
  },
  {
    id: 'multi-round',
    name: '同一 t 多轮：轮次级窗口位置',
    description: '两个窗口从浮动位置开始，每个阶段对象只消耗 1 个 entry；同一个 t 内连续形成 4 轮联合调度并回到基准位置。',
    defaultRbuf: 32,
    times: 4,
    initialPositions: floatedPositions,
    profile: () => ordinary(1, '轻负载普通对象'),
  },
  {
    id: 'independent-early-stop',
    name: '两级独立 FullEarlyStop 追赶',
    description: '两级分别清理自己的全 EarlyStop 链，快路径不消耗普通 entry；遇到第一个普通对象后再进入共享调度。',
    defaultRbuf: 32,
    times: 5,
    initialPositions: floatedPositions,
    profile: ({level, offset}) => {
      if (level === 5 && offset >= 0 && offset <= 1) return fullEarlyStop('Level 5 FullEarlyStop 链');
      if (level === 6 && offset === 0) return fullEarlyStop('Level 6 FullEarlyStop');
      return ordinary(level === 5 ? 2 : 3, '快路径后的普通对象');
    },
  },
];

function findScenario(id) {
  return SCENARIOS.find((scenario) => scenario.id === id) || SCENARIOS[0];
}

function makeGroups(level, demand) {
  const bounded = Math.max(1, Math.min(16, Math.floor(demand)));
  const quotient = Math.floor(bounded / 8);
  const remainder = bounded % 8;
  const groups = [];
  for (let local = 0; local < 8; local += 1) {
    const remaining = quotient + (local < remainder ? 1 : 0);
    if (remaining > 0) {
      groups.push({
        local,
        group: level === 5 ? local : local + 8,
        remaining,
        initial: remaining,
      });
    }
  }
  return groups;
}

function batchText(id) {
  return id >= 0 ? `B${id}` : `B−${Math.abs(id)}`;
}

export class DualWindowSimulator {
  constructor({scenarioId = 'steady', rbuf} = {}) {
    this.scenario = findScenario(scenarioId);
    this.rbuf = evenRows(rbuf ?? this.scenario.defaultRbuf);
    this.rmem = this.rbuf + 2 * WINDOW_ROWS;
    this.memory = [];
    this.stages = new Map();
    this.snapshots = [];
    this.summaries = [];
    this.eventSerial = 0;
    this.nextBatchId = 1;
    this.t = 0;
    this.recovery = {5: 0, 6: 0};
    this.lastIncomingBatch = 0;
    this.lastOutputBatch = null;
    this.forcedHistory = [];
    this.initializeMemory();

    const requested = this.scenario.initialPositions(this.rbuf);
    this.s6 = this.normalizeStart(requested.s6, 0, this.rbuf);
    this.s5 = this.normalizeStart(
      requested.s5,
      this.s6 + WINDOW_ROWS,
      this.rbuf + WINDOW_ROWS,
    );
    if (this.s6 + WINDOW_ROWS > this.s5) this.s5 = this.s6 + WINDOW_ROWS;

    this.roundS5 = this.s5;
    this.roundS6 = this.s6;
    this.initialHeads = {
      5: this.batchAtDecodeRows(5),
      6: this.batchAtDecodeRows(6),
    };
    this.time = null;
  }

  normalizeStart(value, min, max) {
    const even = Math.floor(Number(value) / 2) * 2;
    return Math.max(min, Math.min(max, even));
  }

  initializeMemory() {
    const pairCount = this.rmem / OBJECT_ROWS;
    for (let pair = 0; pair < pairCount; pair += 1) {
      const batch = pair - (pairCount - 1);
      this.memory.push({batch, half: 0}, {batch, half: 1});
    }
  }

  stageKey(level, batch) {
    return `${level}:${batch}`;
  }

  stageProfile(level, batch) {
    const offset = batch - this.initialHeads[level];
    return this.scenario.profile({level, batch, offset, t: this.t});
  }

  getStage(level, batch) {
    const key = this.stageKey(level, batch);
    if (!this.stages.has(key)) {
      const profile = this.stageProfile(level, batch);
      this.stages.set(key, {
        level,
        batch,
        kind: profile.kind,
        note: profile.note,
        groups: profile.kind === 'ordinary' ? makeGroups(level, profile.demand) : [],
        terminal: false,
        terminalReason: null,
        pending: false,
        writebackComplete: false,
        completionT: null,
        serviceCount: 0,
      });
    }
    return this.stages.get(key);
  }

  level5Terminal(batch) {
    if (batch < this.initialHeads[5]) return true;
    const stage = this.stages.get(this.stageKey(5, batch));
    return Boolean(stage && stage.terminal);
  }

  decodeRange(level, startOverride) {
    const start = startOverride ?? (level === 5 ? this.roundS5 : this.roundS6);
    return [start + WINDOW_ROWS - OBJECT_ROWS, start + WINDOW_ROWS];
  }

  batchAtRange([start, end]) {
    if (start < 0 || end > this.memory.length || end - start !== OBJECT_ROWS) return null;
    const first = this.memory[start];
    const second = this.memory[start + 1];
    if (!first || !second || first.batch !== second.batch) return null;
    return first.batch;
  }

  batchAtDecodeRows(level, startOverride) {
    return this.batchAtRange(this.decodeRange(level, startOverride));
  }

  currentStage(level) {
    const batch = this.batchAtDecodeRows(level);
    return batch === null ? null : this.getStage(level, batch);
  }

  canServe(stage) {
    if (!stage || stage.terminal) return false;
    return stage.level !== 6 || this.level5Terminal(stage.batch);
  }

  beginTime() {
    this.roundS5 = this.s5;
    this.roundS6 = this.s6;
    this.time = {
      t: this.t,
      s5Before: this.s5,
      s6Before: this.s6,
      h5Before: this.batchAtDecodeRows(5),
      h6Before: this.batchAtDecodeRows(6),
      inputBatch: this.lastIncomingBatch,
      outputBatch: this.memory[0]?.batch ?? null,
      entries: [],
      rounds: 0,
      c5: 0,
      c6: 0,
      fullEs5: 0,
      fullEs6: 0,
      forced5: [],
      forced6: [],
      events: [],
    };
    this.pushEvent(
      '时刻开始',
      `${batchText(this.lastIncomingBatch)} 已位于 SRAM 底部；从当前 D5/D6 物理位置取得阶段对象。`,
      'start',
    );
  }

  pushEvent(title, detail, phase = 'event', meta = {}) {
    const event = {
      id: this.eventSerial++,
      t: this.t,
      title,
      detail,
      phase,
      ...meta,
    };
    this.time.events.push(event);
    this.recordSnapshot(event);
  }

  recordSnapshot(event) {
    const assertions = this.checkInvariants();
    this.snapshots.push({
      index: this.snapshots.length,
      t: this.t,
      event: clone(event),
      rbuf: this.rbuf,
      rmem: this.rmem,
      memory: clone(this.memory),
      s5: this.roundS5,
      s6: this.roundS6,
      committedS5: this.s5,
      committedS6: this.s6,
      timeStartS5: this.time.s5Before,
      timeStartS6: this.time.s6Before,
      nextS5: Number.isFinite(event.nextS5) ? event.nextS5 : null,
      nextS6: Number.isFinite(event.nextS6) ? event.nextS6 : null,
      d5: this.decodeRange(5),
      d6: this.decodeRange(6),
      h5: this.batchAtDecodeRows(5),
      h6: this.batchAtDecodeRows(6),
      recovery: clone(this.recovery),
      entries: clone(this.time.entries),
      rounds: this.time.rounds,
      remainingEntries: ENTRY_BUDGET - this.time.entries.length,
      c5: this.time.c5,
      c6: this.time.c6,
      fullEs5: this.time.fullEs5,
      fullEs6: this.time.fullEs6,
      forced5: clone(this.time.forced5),
      forced6: clone(this.time.forced6),
      incomingBatch: this.lastIncomingBatch,
      outputBatch: this.memory[0]?.batch ?? null,
      assertions,
      events: clone(this.time.events),
    });
  }

  positionLegal(s5, s6) {
    return s6 >= 0 && s6 + WINDOW_ROWS <= s5 && s5 + WINDOW_ROWS <= this.rmem;
  }

  applyRecovery() {
    const before = {s5: this.roundS5, s6: this.roundS6};
    const moved = [];
    let progress = true;
    while (progress) {
      progress = false;
      const want5 = this.recovery[5] >= OBJECT_ROWS;
      const want6 = this.recovery[6] >= OBJECT_ROWS;

      if (want5 && want6 && this.positionLegal(this.roundS5 + 2, this.roundS6 + 2)) {
        this.roundS5 += 2;
        this.roundS6 += 2;
        this.recovery[5] -= 2;
        this.recovery[6] -= 2;
        moved.push('L5', 'L6');
        progress = true;
        continue;
      }
      if (want5 && this.positionLegal(this.roundS5 + 2, this.roundS6)) {
        this.roundS5 += 2;
        this.recovery[5] -= 2;
        moved.push('L5');
        progress = true;
        continue;
      }
      if (want6 && this.positionLegal(this.roundS5, this.roundS6 + 2)) {
        this.roundS6 += 2;
        this.recovery[6] -= 2;
        moved.push('L6');
        progress = true;
      }
    }
    return {
      before,
      after: {s5: this.roundS5, s6: this.roundS6},
      moved,
      blocked5: this.recovery[5] > 0,
      blocked6: this.recovery[6] > 0,
    };
  }

  completeStage(stage, reason) {
    if (stage.terminal) throw new Error(`stage completed twice: ${stage.level}:${stage.batch}`);
    stage.terminal = true;
    stage.terminalReason = reason;
    stage.pending = false;
    stage.writebackComplete = true;
    stage.completionT = this.t;
    this.recovery[stage.level] += OBJECT_ROWS;
    if (stage.level === 5) {
      this.time.c5 += 1;
      if (reason === 'FullEarlyStop') this.time.fullEs5 += 1;
    } else {
      this.time.c6 += 1;
      if (reason === 'FullEarlyStop') this.time.fullEs6 += 1;
    }
  }

  clearFullEarlyStop() {
    let any = false;
    let guard = 0;
    while (guard++ < this.rmem) {
      const candidates = [];
      for (const level of [5, 6]) {
        if (this.recovery[level] > 0) continue;
        const stage = this.currentStage(level);
        if (stage && this.canServe(stage) && stage.kind === 'full-es') candidates.push(stage);
      }
      if (candidates.length === 0) break;
      any = true;
      for (const stage of candidates) this.completeStage(stage, 'FullEarlyStop');
      const movement = this.applyRecovery();
      const names = candidates.map((stage) => `L${stage.level}(${batchText(stage.batch)})`).join('、');
      const blocked = [
        movement.blocked5 ? 'L5 推进受阻' : '',
        movement.blocked6 ? 'L6 推进受阻' : '',
      ].filter(Boolean).join('；');
      this.pushEvent(
        'FullEarlyStop 清理',
        `${names} 执行原 EarlyStopAction 并结束阶段，不消耗普通 entry。${blocked ? ` ${blocked}，推进欠账保留。` : ''}`,
        'early-stop',
        {levels: candidates.map((stage) => stage.level)},
      );
      if (candidates.every((stage) => this.recovery[stage.level] > 0)) break;
    }
    return any;
  }

  ordinaryCandidates() {
    const candidates = [];
    for (const level of [5, 6]) {
      if (this.recovery[level] > 0) continue;
      const stage = this.currentStage(level);
      if (!stage || !this.canServe(stage) || stage.kind !== 'ordinary' || stage.terminal) continue;
      for (const group of stage.groups) {
        if (group.remaining > 0) candidates.push({stage, group});
      }
    }
    candidates.sort((a, b) => {
      if (b.group.remaining !== a.group.remaining) return b.group.remaining - a.group.remaining;
      if (a.stage.level !== b.stage.level) return a.stage.level - b.stage.level;
      return a.group.group - b.group.group;
    });
    return candidates;
  }

  executeOrdinaryRound() {
    const remainingBudget = ENTRY_BUDGET - this.time.entries.length;
    if (remainingBudget <= 0) return false;
    const candidates = this.ordinaryCandidates();
    if (candidates.length === 0) return false;

    const selected = candidates.slice(0, remainingBudget);
    this.time.rounds += 1;
    const round = this.time.rounds;
    const touched = new Set();
    for (const candidate of selected) {
      const {stage, group} = candidate;
      const visit = group.initial - group.remaining;
      group.remaining -= 1;
      stage.serviceCount += 1;
      stage.pending = true;
      touched.add(stage);
      const codeBase = group.group * 4;
      const hiso = codeBase + (visit % 2);
      const siso = codeBase + 2 + (visit % 2);
      this.time.entries.push({
        slot: this.time.entries.length,
        round,
        level: stage.level,
        batch: stage.batch,
        group: group.group,
        hiso,
        siso,
      });
    }

    const completed = [];
    for (const stage of touched) {
      if (stage.groups.every((group) => group.remaining === 0)) {
        this.completeStage(stage, 'Normal');
        completed.push(stage);
      }
    }
    const movement = this.applyRecovery();
    const allocation = selected
      .map(({stage, group}) => `E${this.time.entries.find((entry) => entry.round === round && entry.level === stage.level && entry.group === group.group)?.slot ?? '?'}:L${stage.level}/${batchText(stage.batch)}/G${group.group}`)
      .join('，');
    const completionText = completed.length > 0
      ? ` 完成：${completed.map((stage) => `L${stage.level}(${batchText(stage.batch)})`).join('、')}。`
      : ' 当前对象仍有 pending。';
    const debtText = this.recovery[5] || this.recovery[6]
      ? ` 推进欠账：L5=${this.recovery[5]} 行，L6=${this.recovery[6]} 行。`
      : '';
    this.pushEvent(
      `普通联合调度 Round ${round}`,
      `${allocation}。完成本轮执行和 SRAM 写回。${completionText}${debtText}`,
      'ordinary',
      {round, completed: completed.map((stage) => stage.level)},
    );
    return true;
  }

  forceStage(level, reason, stageOverride = null) {
    const stage = stageOverride || this.currentStage(level);
    if (!stage || stage.terminal) return null;
    stage.terminal = true;
    stage.terminalReason = 'ForcedEvicted';
    stage.pending = false;
    stage.writebackComplete = false;
    stage.completionT = this.t;
    const record = {level, batch: stage.batch, t: this.t, reason};
    this.forcedHistory.push(record);
    this.time[`forced${level}`].push(stage.batch);
    return record;
  }

  resolveBoundary() {
    const pre = {s5: this.roundS5, s6: this.roundS6};
    const current5 = this.currentStage(5);
    const current6 = this.currentStage(6);
    let next5 = this.roundS5 - OBJECT_ROWS;
    let next6 = this.roundS6 - OBJECT_ROWS;
    const forced = [];

    if (next6 < 0) {
      if (this.recovery[6] >= OBJECT_ROWS) {
        this.recovery[6] -= OBJECT_ROWS;
        next6 = this.roundS6;
      } else {
        const record = this.forceStage(
          6,
          'Level 6 到达 SRAM 顶部，无法继续跟随 pending 对象',
          current6,
        );
        if (record) forced.push(record);
        next6 = 0;
      }
    }

    if (next5 < 0) next5 = 0;
    this.roundS5 = next5;
    this.roundS6 = next6;

    // 必要跟随位置先确定；完成对象产生的恢复只使用剩余合法空间。
    this.applyRecovery();

    if (this.roundS6 + WINDOW_ROWS > this.roundS5) {
      if (current5 && !current5.terminal) {
        const record = this.forceStage(
          5,
          'Level 5 必要跟随在联合裁决后仍与 Level 6 重叠',
          current5,
        );
        if (record) forced.push(record);
      }
      this.recovery[5] = 0;
      this.roundS5 = this.roundS6 + WINDOW_ROWS;
    }

    if (this.roundS5 + WINDOW_ROWS > this.rmem) {
      this.roundS5 = this.rmem - WINDOW_ROWS;
    }

    const forcedText = forced.length > 0
      ? ` 强制结束：${forced.map((item) => `L${item.level}(${batchText(item.batch)})`).join('、')}；物理数据仍保留在 SRAM。`
      : ' 未触发 forced-evicted。';
    const resolved = {s5: this.roundS5, s6: this.roundS6};
    // 边界事件仍展示上推前的 SRAM 和窗口；resolved 位置只在物理上推后提交。
    this.roundS5 = pre.s5;
    this.roundS6 = pre.s6;
    this.pushEvent(
      '时刻边界联合裁决',
      `必要跟随优先，恢复只使用剩余合法空间。下一位置 L6 ${pre.s6}→${resolved.s6}，L5 ${pre.s5}→${resolved.s5}。${forcedText}`,
      'boundary',
      {forced: clone(forced), nextS5: resolved.s5, nextS6: resolved.s6},
    );

    const output = this.memory.slice(0, OBJECT_ROWS);
    const outputBatch = output[0]?.batch ?? null;
    this.memory = this.memory.slice(OBJECT_ROWS);
    const incomingBatch = this.nextBatchId++;
    this.memory.push({batch: incomingBatch, half: 0}, {batch: incomingBatch, half: 1});

    this.s5 = resolved.s5;
    this.s6 = resolved.s6;
    this.lastOutputBatch = outputBatch;
    this.lastIncomingBatch = incomingBatch;

    const summary = {
      t: this.t,
      inputBatch: this.time.inputBatch,
      outputBatch,
      h5Before: this.time.h5Before,
      h6Before: this.time.h6Before,
      s5Before: this.time.s5Before,
      s6Before: this.time.s6Before,
      s5After: this.s5,
      s6After: this.s6,
      c5: this.time.c5,
      c6: this.time.c6,
      fullEs5: this.time.fullEs5,
      fullEs6: this.time.fullEs6,
      forced5: clone(this.time.forced5),
      forced6: clone(this.time.forced6),
      entries: this.time.entries.length,
      rounds: this.time.rounds,
      recovery5: this.recovery[5],
      recovery6: this.recovery[6],
      finalSnapshotIndex: this.snapshots.length - 1,
    };
    this.summaries.push(summary);

    this.t += 1;
    this.roundS5 = this.s5;
    this.roundS6 = this.s6;
    this.time = {
      t: this.t,
      s5Before: this.s5,
      s6Before: this.s6,
      h5Before: this.batchAtDecodeRows(5),
      h6Before: this.batchAtDecodeRows(6),
      inputBatch: incomingBatch,
      outputBatch,
      entries: [],
      rounds: 0,
      c5: 0,
      c6: 0,
      fullEs5: 0,
      fullEs6: 0,
      forced5: [],
      forced6: [],
      events: [],
    };
    this.pushEvent(
      '物理推进完成',
      `${batchText(outputBatch)} 的两行从顶部输出；SRAM 整体上移 2 行；${batchText(incomingBatch)} 写入底部。`,
      'shift',
      {fromT: this.t - 1, incomingBatch, outputBatch},
    );
  }

  runTime() {
    if (!this.time) this.beginTime();
    this.clearFullEarlyStop();
    let guard = 0;
    while (this.time.entries.length < ENTRY_BUDGET && guard++ < ENTRY_BUDGET + 2) {
      const progressed = this.executeOrdinaryRound();
      if (!progressed) break;
      this.clearFullEarlyStop();
    }
    // 普通预算耗尽以后，仍允许清理本轮写回后新暴露的 FullEarlyStop 链。
    this.clearFullEarlyStop();
    this.resolveBoundary();
  }

  simulate(times = this.scenario.times) {
    const count = Math.max(1, Math.min(20, Math.floor(times)));
    while (this.summaries.length < count) this.runTime();
    return {
      snapshots: this.snapshots,
      summaries: this.summaries,
      scenario: this.scenario,
      constants: {WINDOW_ROWS, OBJECT_ROWS, ENTRY_BUDGET, BLOCK_COLUMNS},
    };
  }

  checkInvariants() {
    const checks = [];
    const add = (name, pass, detail) => checks.push({name, pass: Boolean(pass), detail});
    add('窗口位于 SRAM 内', this.roundS6 >= 0 && this.roundS5 + WINDOW_ROWS <= this.rmem,
      `0 <= ${this.roundS6}；${this.roundS5}+22 <= ${this.rmem}`);
    add('两个窗口不重叠', this.roundS6 + WINDOW_ROWS <= this.roundS5,
      `${this.roundS6}+22 <= ${this.roundS5}`);
    add('窗口保持 2 行对齐', this.roundS5 % 2 === 0 && this.roundS6 % 2 === 0,
      `S5=${this.roundS5}，S6=${this.roundS6}`);
    add('D5 映射到一个完整 B', this.batchAtDecodeRows(5) !== null,
      `D5=[${this.decodeRange(5).join(',')})`);
    add('D6 映射到一个完整 B', this.batchAtDecodeRows(6) !== null,
      `D6=[${this.decodeRange(6).join(',')})`);
    const used = this.time?.entries.length ?? 0;
    add('普通 entry 不超过 8', used <= ENTRY_BUDGET, `used=${used}`);
    const badL6 = (this.time?.entries ?? []).find((entry) => entry.level === 6 && !this.level5Terminal(entry.batch));
    add('Level 6 不早于 Level 5', !badL6,
      badL6 ? `${batchText(badL6.batch)} 的 L5 尚未结束` : '所有已服务 L6 对象均满足先后关系');
    return checks;
  }
}

export function runSelfTests() {
  const results = [];
  const verify = (name, fn) => {
    try {
      fn();
      results.push({name, pass: true});
    } catch (error) {
      results.push({name, pass: false, detail: error.message});
    }
  };
  const assert = (condition, message) => {
    if (!condition) throw new Error(message);
  };

  for (const scenario of SCENARIOS) {
    verify(`${scenario.name}：全部不变量`, () => {
      const sim = new DualWindowSimulator({scenarioId: scenario.id, rbuf: scenario.defaultRbuf});
      const trace = sim.simulate(Math.min(4, scenario.times));
      const failure = trace.snapshots.flatMap((snapshot) => snapshot.assertions).find((item) => !item.pass);
      assert(!failure, failure ? `${failure.name}: ${failure.detail}` : 'unknown');
    });
  }

  verify('双 pending 相邻窗口同步移动', () => {
    const sim = new DualWindowSimulator({scenarioId: 'both-pending', rbuf: 32});
    const {summaries} = sim.simulate(1);
    assert(summaries[0].s5After === summaries[0].s5Before - 2, 'Level 5 未向顶部移动 2 行');
    assert(summaries[0].s6After === summaries[0].s6Before - 2, 'Level 6 未向顶部移动 2 行');
    assert(summaries[0].forced5.length === 0, 'Level 5 被错误强制弹出');
  });

  verify('多轮场景使用轮次级位置', () => {
    const sim = new DualWindowSimulator({scenarioId: 'multi-round', rbuf: 32});
    const {summaries} = sim.simulate(1);
    assert(summaries[0].rounds === 4, `期望 4 轮，实际 ${summaries[0].rounds}`);
    assert(summaries[0].entries === 8, `期望使用 8 entry，实际 ${summaries[0].entries}`);
    assert(summaries[0].c5 === 4 && summaries[0].c6 === 4, '两级完成量不是 4/4');
  });

  verify('Level 6 forced-evicted 不删除物理数据', () => {
    const sim = new DualWindowSimulator({scenarioId: 'top-boundary', rbuf: 32});
    const forcedBatch = sim.initialHeads[6];
    const {summaries, snapshots} = sim.simulate(1);
    assert(summaries[0].forced6.includes(forcedBatch), '没有触发预期的 Level 6 forced-evicted');
    const finalMemory = snapshots[snapshots.length - 1].memory;
    assert(finalMemory.some((row) => row.batch === forcedBatch), 'forced-evicted 数据被错误删除');
  });

  return results;
}

export const SIM_CONSTANTS = {WINDOW_ROWS, OBJECT_ROWS, ENTRY_BUDGET, BLOCK_COLUMNS};
