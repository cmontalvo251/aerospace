function write_to_csfile(freq)

newupdaterate = [num2str(freq),'	!Feedback Update Rate(Hz)'];

data = {'1	!Aircraft Control on or off';
'-0.02	!Kpy Proportional gain on crossrange (aircraft)';
'0.1	!Kdy Derivative gain on crossrange (aircraft)';
'-0.8	!Kpsi Proportional gain on heading (aircraft)';
'1.2	!Kpp proportional gain on roll angle (aircraft)';
'0.1	!Kdp Derivative gain on roll angle (aircraft)';
'-0.06	!Kpz Proportional gain on altitude (aircraft)';
'0.02	!Kdz Derivative gain on altitude (aircraft)';
'-1.0	!Kpt Proportional gain on pitch (aircraft)';
'-0.15	!Kv proportional gain on sideslip (aircraft)';
'2.0	!Ku proportional gain on speed (aircraft)    ';
newupdaterate;
'1	!Turn On Sensor Errors'};

writedata('Input_Files/SBXC.CS',data);