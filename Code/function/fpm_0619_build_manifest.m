function manifest = fpm_0619_build_manifest(cfg, stage, extra)
%FPM_0619_BUILD_MANIFEST 生成统一的追踪信息 schema
% 用法：
%   manifest = fpm_0619_build_manifest(cfg, "preproc", extra_struct);
%   manifest = fpm_0619_build_manifest(cfg, "recover", extra_struct);
%
% 统一字段设计：
%   manifest.schema_version
%   manifest.project
%   manifest.batch_folder
%   manifest.stage
%   manifest.preset_name
%   manifest.source
%   manifest.optics
%   manifest.calibration
%   manifest.geometry
%   manifest.recon
%   manifest.outputs
%   manifest.notes
%
% 这样预处理和重建都用同一套顶层结构，减少追踪时的分叉。

if nargin < 3 || isempty(extra)
    extra = struct();
end

manifest = struct();
manifest.schema_version = "fpm_manifest_v1";
manifest.project = "FPM_0619";
manifest.batch_folder = cfg.batch_folder;
manifest.stage = string(stage);
manifest.preset_name = cfg.preset_name;

manifest.source = struct();
if isfield(cfg, 'input_subdir')
    manifest.source.input_subdir = cfg.input_subdir;
else
    manifest.source.input_subdir = "";
end
if isfield(cfg, 'data_name')
    manifest.source.data_name = cfg.data_name;
else
    manifest.source.data_name = "";
end
if isfield(cfg, 'preproc_output_tag')
    manifest.source.preproc_output_tag = cfg.preproc_output_tag;
else
    manifest.source.preproc_output_tag = "";
end
if isfield(cfg, 'preproc_ptr')
    manifest.source.preproc_ptr = cfg.preproc_ptr;
else
    manifest.source.preproc_ptr = "";
end

manifest.optics = cfg.optics;
manifest.calibration = cfg.calibration;
manifest.geometry = cfg.geometry;
manifest.recon = cfg.recon;
manifest.outputs = struct();
manifest.notes = struct();

extra_fields = fieldnames(extra);
for k = 1:numel(extra_fields)
    manifest.(extra_fields{k}) = extra.(extra_fields{k});
end
end
