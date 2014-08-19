function cfg = run_ft_artifact(cfg, data, typename)
%   RUN_FT_ARTIFACT   Runs fieldtrip artifact rejection
%     [CFG] = RUN_FT_ARTIFACT(CFG, DATA, TYPENAME)
% 
%   TODO
%   
%   Created by Zarrar Shehzad on 2012-09-11.
%

if isfield(cfg.artfctdef.(typename), 'artifact')
    prev_artifact = cfg.artfctdef.(typename).artifact;
    cfg.artfctdef.(typename) = rmfield(cfg.artfctdef.(typename), 'artifact');
else
    prev_artifact = [];
end

[cfg cur_artifact] = eval(['ft_artifact_' typename, '(cfg, data)']);

if ~isempty(prev_artifact)
    inds = arrayfun(@(x) any(x == prev_artifact(:,1)), cur_artifact(:,1));
    cfg.artfctdef.(typename).artifact = [prev_artifact; cur_artifact(inds,:)];
end

end %  function
