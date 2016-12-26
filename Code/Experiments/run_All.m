function run_All()


matchingMethod=1;
tran=1000;
reps=10;
run_Swiss_Roll(matchingMethod,tran,reps);
tran=500;
run_Swiss_Roll_Noise(matchingMethod,tran,reps);
run_Swiss_Roll_Outlier(matchingMethod,tran,reps);