function run_golden_master(mode)
% run_golden_master  Регрессионная обвязка Golden Master для ядра LoRaPHYmod.
%
% Назначение:
%   Фиксирует эталонное поведение вычислительного ядра (BER/PER полного тракта
%   и метрики ML-оценки канала) при фиксированном seed и проверяет его
%   неизменность после рефакторинга. Служит защитным барьером (merge gate):
%   любое изменение числового результата обнаруживается до слияния ветки.
%
% Режимы:
%   run_golden_master record   — исполнить кейсы и записать эталоны в tests/golden/.
%                                Запускать ОДИН раз на опорном коммите
%                                (baseline-pre-refactor), до любых правок ядра.
%   run_golden_master verify   — исполнить кейсы и сверить с эталонами.
%                                Запускать после каждого шага рефакторинга.
%                                При расхождении завершается ошибкой.
%
% По умолчанию (без аргумента) — режим verify.
%
% Гранулярность:
%   Обвязка гарантирует КОМПОНЕНТЫ ядра (модем, каналы, симулятор,
%   estimateChannel), а не тела сценариев. Кейсы заданы как данные —
%   эскиз паттерна spec-as-data для будущего обобщённого sweep-драйвера.
%
% Строгость сверки:
%   Кейсы 'link'  : tol = 0  — точное совпадение (счётчики целочисленны,
%                             BER/PER детерминированы).
%   Кейсы 'chanest': tol = 1e-12 (отн.) — плавающие агрегаты; запас на
%                    возможные межплатформенные различия суммирования.
%
% Расположение эталонов:
%   tests/golden/<name>.mat  (вне results/ ⇒ отслеживается git и коммитится).
%
% Зависимости: gm_run_case.m (в той же папке tests/).
% Совместимость: MATLAB R2024a.
% ---------------------------------------------------------------

    if nargin < 1 || isempty(mode)
        mode = 'verify';
    end
    mode = lower(string(mode));

    tests_dir   = fileparts(mfilename('fullpath'));
    projectRoot = fileparts(tests_dir);
    addpath(genpath(projectRoot));

    golden_dir = fullfile(tests_dir, 'golden');
    if ~exist(golden_dir, 'dir')
        mkdir(golden_dir);
    end

    cases = gm_define_cases();

    switch mode
        case "record"
            gm_do_record(cases, golden_dir, projectRoot);
        case "verify"
            gm_do_verify(cases, golden_dir);
        otherwise
            error('run_golden_master:badMode', ...
                'Режим должен быть ''record'' или ''verify'' (получено: %s).', mode);
    end
end

% ------------------------------------------------------------------
%  Определения кейсов (spec-as-data)
% ------------------------------------------------------------------
function cases = gm_define_cases()
    % Общие PHY-параметры (EBYTE E22-868T22U, SX1262, EU868)
    common = struct( ...
        'rf_freq',     868e6, ...
        'sf',          7, ...
        'bw',          125e3, ...
        'fs',          1e6, ...
        'payloadBits', 128);

    % --- Кейс 1: полный тракт, канал Релея TDL, режим детектирования (CR=1) ---
    % Гарантирует связку модем + RayleighTDLChannel (PerSymbol/AR(1)) + симулятор.
    c1 = common;
    c1.name        = "link_tdl";
    c1.kind        = 'link';
    c1.channel     = 'tdl';
    c1.seed        = 20260711;
    c1.tol         = 0;
    c1.CR          = 1;
    c1.PreambleLen = 8;
    c1.FastMode    = true;
    c1.Npkts       = 24;
    c1.snr_list    = [-8, -4, 0, 4];
    c1.v_mps_list  = [0, 30];
    c1.tdlDelays   = [0, 0.5e-6, 2.0e-6, 4.0e-6];
    c1.tdlGains    = [0, -4,     -10,    -18];

    % --- Кейс 2: полный тракт, АБГШ, режим коррекции (CR=4, Хэмминг) ---
    % Гарантирует ветвь исправляющего декодера, отличную от кейса 1.
    c2 = common;
    c2.name        = "link_awgn";
    c2.kind        = 'link';
    c2.channel     = 'awgn';
    c2.seed        = 20260712;
    c2.tol         = 0;
    c2.CR          = 4;
    c2.PreambleLen = 8;
    c2.FastMode    = true;
    c2.Npkts       = 24;
    c2.snr_list    = [-10, -8, -6, -4];
    c2.v_mps_list  = 0;                 % для АБГШ скорость не влияет
    c2.tdlDelays   = [];                % не используются
    c2.tdlGains    = [];

    % --- Кейс 3: ML/MMSE-оценка канала (estimateChannel) ---
    % Гарантирует числовой выход оценщика (RMSE амплитуды, дисперсия).
    c3 = common;
    c3.name             = "chanest";
    c3.kind             = 'chanest';
    c3.seed             = 20260713;
    c3.tol              = 1e-12;
    c3.CR               = 1;
    c3.snr_list         = [-6, 0, 6];
    c3.preamble_configs = [1, 8];
    c3.N_mc             = 200;

    cases = {c1, c2, c3};
end

% ------------------------------------------------------------------
%  Режим record — запись эталонов
% ------------------------------------------------------------------
function gm_do_record(cases, golden_dir, projectRoot)
    fprintf('\n=== Golden Master: ЗАПИСЬ эталонов ===\n');
    meta = gm_meta(projectRoot);
    fprintf('MATLAB %s | git %s | %s\n\n', ...
        meta.matlab_version, meta.git_hash, meta.timestamp);

    for k = 1:numel(cases)
        cs  = cases{k};
        fprintf('  [%d/%d] %-12s ... ', k, numel(cases), cs.name);
        t0  = tic;
        res = gm_run_case(cs);            % #ok<NASGU>  сохраняется в ref
        ref = struct('name', cs.name, 'kind', cs.kind, 'seed', cs.seed, ...
                     'tol', cs.tol, 'spec', cs, 'data', res, 'meta', meta);
        fname = fullfile(golden_dir, char(cs.name) + ".mat");
        save(fname, 'ref');
        fprintf('готово (%.1f с) → %s\n', toc(t0), fname);
    end

    fprintf('\nЭталоны записаны. Не забудьте закоммитить tests/golden/*.mat.\n\n');
end

% ------------------------------------------------------------------
%  Режим verify — сверка с эталонами
% ------------------------------------------------------------------
function gm_do_verify(cases, golden_dir)
    fprintf('\n=== Golden Master: СВЕРКА с эталонами ===\n\n');
    fprintf('%-12s  %-14s  %-12s  %-6s\n', 'Кейс', 'Поле', 'max|Δ|', 'Статус');
    fprintf('%s\n', repmat('-', 1, 52));

    n_fail = 0;

    for k = 1:numel(cases)
        cs    = cases{k};
        fname = fullfile(golden_dir, char(cs.name) + ".mat");

        if ~exist(fname, 'file')
            fprintf('%-12s  ЭТАЛОН ОТСУТСТВУЕТ (%s)\n', cs.name, fname);
            n_fail = n_fail + 1;
            continue;
        end

        S   = load(fname, 'ref');
        ref = S.ref;
        cur = gm_run_case(cs);

        flds = fieldnames(ref.data);
        for fi = 1:numel(flds)
            f       = flds{fi};
            a       = ref.data.(f);
            b       = cur.(f);
            [ok, dev] = gm_compare(a, b, cs.tol);
            status  = "OK";
            if ~ok
                status = "FAIL";
                n_fail = n_fail + 1;
            end
            fprintf('%-12s  %-14s  %-12.3e  %-6s\n', cs.name, f, dev, status);
        end
    end

    fprintf('%s\n', repmat('-', 1, 52));
    if n_fail == 0
        fprintf('РЕЗУЛЬТАТ: все кейсы совпали с эталоном (GREEN).\n\n');
    else
        error('run_golden_master:regression', ...
            'РЕГРЕССИЯ: расхождений — %d. Слияние блокируется.', n_fail);
    end
end

% ------------------------------------------------------------------
%  Сверка двух числовых массивов
% ------------------------------------------------------------------
function [ok, dev] = gm_compare(a, b, tol)
    % Возвращает признак совпадения ok и максимальное отклонение dev.
    if ~isequal(size(a), size(b))
        ok = false; dev = Inf; return;
    end
    dev = max(abs(a(:) - b(:)), [], 'omitnan');
    if isempty(dev), dev = 0; end

    if tol == 0
        ok = isequaln(a, b);            % бит-в-бит (с корректной обработкой NaN)
    else
        scale = max(1, max(abs(a(:)), [], 'omitnan'));
        ok    = dev <= tol * scale;     % относительный допуск
    end
end

% ------------------------------------------------------------------
%  Метаданные провенанса
% ------------------------------------------------------------------
function meta = gm_meta(projectRoot)
    meta.matlab_version = version;
    meta.timestamp      = datestr(now, 'yyyy-mm-dd HH:MM:SS'); %#ok<TNOW1,DATST>

    cmd = sprintf('git -C "%s" rev-parse --short HEAD', projectRoot);
    [st, out] = system(cmd);
    if st == 0
        meta.git_hash = strtrim(out);
    else
        meta.git_hash = '(git недоступен)';
    end
end
