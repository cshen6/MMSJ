function run_All(matchingMethod)

if nargin<1
matchingMethod=1;
end
tran=1000;
reps=20;
run_Swiss_Roll(1,matchingMethod,tran,reps);
run_Swiss_Roll(2,matchingMethod,tran,reps);
run_Swiss_Roll(3,matchingMethod,tran,reps);