function data = read_quantized_llr(paths, makePlot)
%READ_QUANTIZED_LLR  读取 C++ 导出的三类 LLR 并可选画直方图。
%   DATA = READ_QUANTIZED_LLR()  从默认路径读取三份文件：
%       data/quantized_llr_float.txt        (量化前)
%       data/quantized_llr_quantized_codes.txt (码值)
%       data/quantized_llr_dequantized.txt   (反量化后)
%   DATA = READ_QUANTIZED_LLR(PATHS)  PATHS 为 struct，字段：
%       .float, .codes, .dequantized  （任意缺失则用默认路径）
%   DATA = READ_QUANTIZED_LLR(..., MAKEPLOT) 当 MAKEPLOT 为 true 时绘制 3x1 子图直方图。
%
%   返回 DATA 结构体：DATA.float / DATA.codes / DATA.dequantized。

if nargin < 1 || isempty(paths)
  paths = struct();
end
if nargin < 2
  makePlot = false;
end

thisDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(thisDir);  % ../../

% 默认路径
default.float       = fullfile(repoRoot, 'data', 'llr', 'quantized_llr_float.txt');
default.codes       = fullfile(repoRoot, 'data', 'llr', 'quantized_llr_quantized_codes.txt');
default.dequantized = fullfile(repoRoot, 'data', 'llr', 'quantized_llr_dequantized.txt');

% 合并路径
fields = {'float','codes','dequantized'};
for i = 1:numel(fields)
  key = fields{i};
  if ~isfield(paths, key) || isempty(paths.(key))
    paths.(key) = default.(key);
  end
end

% 读取函数
readFile = @(p,label,suffix) local_read(p, label, suffix);

data.float       = readFile(paths.float, 'float', 'float');
data.codes       = readFile(paths.codes, 'codes', 'quantized_codes');
data.dequantized = readFile(paths.dequantized, 'dequantized', 'dequantized');

if makePlot
  figure;
  tiledlayout(3,1, 'TileSpacing','compact');

  nexttile;
  histogram(data.float, 'Normalization','pdf');
  xlabel('LLR'); ylabel('PDF'); title(sprintf('LLR (float, n=%d)', numel(data.float)));

  nexttile;
  histogram(data.codes, 'Normalization','pdf');
  xlabel('Quantized code'); ylabel('PDF'); title(sprintf('Quantized codes (n=%d)', numel(data.codes)));

  nexttile;
  histogram(data.dequantized, 'Normalization','pdf');
  xlabel('LLR'); ylabel('PDF'); title(sprintf('LLR (dequantized, n=%d)', numel(data.dequantized)));
end
end

function v = local_read(path, label, suffix)
if ~isfile(path)
  % 尝试按照默认前缀/后缀在同目录下匹配 quantized_llr_*_suffix.txt
  [dirp,~,ext] = fileparts(path);
  pattern = fullfile(dirp, sprintf('quantized_llr_*_%s%s', suffix, ext));
  files = dir(pattern);
  if ~isempty(files)
    path = fullfile(files(1).folder, files(1).name);
  else
    error('read_quantized_llr:FileNotFound', '找不到 %s 文件：%s', label, path);
  }
end
v = readmatrix(path, 'FileType','text');
end
