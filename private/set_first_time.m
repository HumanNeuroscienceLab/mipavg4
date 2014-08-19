function data = set_first_time(data, new_first_pt, logg)
%   SET_FIRST_TIME   Change the first point in time in SECONDS.
%     [DATA] = SET_FIRST_TIME(DATA, new_first_pt, logg)
% 
%   This function uses ft_redefinetrial and its cfg.offset option
%   to change the time-axis.
%   
%   Created by Zarrar Shehzad on 2012-09-24.

    old_first_pt  = data.time{1}(1);
    time_diff     = new_first_pt - old_first_pt;
    
    if time_diff
        logg.write('\nAdjusting start of time-axis from %.4f to %.4f secs\n',  ...
                    old_first_pt, new_first_pt);
        
        if mod(new_first_pt, 1/data.fsample)
            logg.write(['Note: new time point isn''t a multiple of the ' ...
                        'sampling rate. There will be some rounding.\n'])
        end
        
        tmpevent    = data.cfg.event;
        tmptrl      = data.cfg.trl;
        tmptrlold   = data.cfg.trlold;
                    
        cfg = [];
        cfg.offset = round(time_diff * data.fsample);        
        data = ft_redefinetrial(cfg, data);
        
        data.cfg.event  = tmpevent;
        data.cfg.trl    = tmptrl;
        data.cfg.trlold = tmptrlold;
        clear tmpevent tmptrl tmptrlold;
        
        logg.write('\nFirst time point is now %.4f secs\n', data.time{1}(1));
    else
        logg.write('\nNOT adjusting start of time-axis %.4f secs',  ...
                    old_first_pt);
    end
    
end %  function