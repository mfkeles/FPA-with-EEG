function M = define_M()
valueSet = {'REM','NREM','Active_Wake','Quiet_Wake'};
keySet = [3,2,4,1];
M =  containers.Map(keySet,valueSet);
end