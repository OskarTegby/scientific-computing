%-- Play Movie

if (size(sol,1)==0)
  disp(['Run solver!',char(7)]);
  return;
end
if (size(fmat,1)==0)
  disp(['Record Movie!',char(7)]);
  return;
end

nrepeat = str2num(get(repeat_ui,'String'));
pics = get(speed_ui,'Value');
movie(f1,fmat,nrepeat-1,pics);