MRJ: For the moment, adding this file simply to verify that my GitHub access is working correctly.

But every project should have a README, right?

Perhaps in the future, this document will have more useful contents such as a version history, usage notes, other documentation, or 
 other info?

To start off, maybe I'll list MRJ's aims for commits in the near future:
* Fix bug (?) where including/excluding certain stimulus codes in SGC files can yield incorrect #s of occurrences for OTHER stimulus 
  codes (which in theory should have nothing to do with the codes included/excluded).
* * That problem is now fixed -- but still need to track down another unexplained phenomenon wherein including/excluding certain 
    condition codes from analysis affects the numbers of OTHER conditons that survive artifact rejection.
* Add in ability to linearly detrend epochs (which can be approximated right now with the right high-pass filter, I guess, but that 
  isn't quite the same...)
* Add in ability to regress out HEOG/VEOG (or potentially an arbitrary set of other channels) from all other channels' data to 
  remove any residual influence from small eye movements.

# Known Issues

## supergui

On some machines, you might see this error: "Error using supergui (line 126) supergui error: argument 'fig' must be numeric".

According to this [website](http://sccn.ucsd.edu/pipermail/eeglablist/2014/008851.html) you just need to change line 126 in supergui from

'fig'       'real'       []      0;

to 

'fig'       ''       []      0;
