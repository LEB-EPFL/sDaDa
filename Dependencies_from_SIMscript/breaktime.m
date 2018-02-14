function [weeks, days, hours, mins, secs]=breaktime(totseconds, varargin)
%breaktime breaks a number of seconds into hours, minutes and seconds
%Format: [weeks, days, hours, mins, secs]= breaktime(totalseconds)

weeks=floor(totseconds/604800);
remsecs=rem(totseconds,604800);
days=floor(remsecs/86400);
remsecs=rem(remsecs,86400);
hours=floor(remsecs/3600);
remsecs=rem(remsecs,3600);
mins=floor(remsecs/60);
secs=rem(remsecs,60);
if ~isempty('varargin')
    if strcmp(varargin, 'disp')
        disp([num2str(weeks) ' weeks, ' num2str(days) ' days, ' num2str(hours) ' hours, ' num2str(mins) ' minutes and ' num2str(secs) ' seconds.'])
    end
end
end