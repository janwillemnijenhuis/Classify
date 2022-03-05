// SPECIFY TO THE LOCATION OF YOUR FILES // 
cd "C:\Users\janwi\OneDrive\Documents\Classify\Classify\Old Swopit commands" 
mata: mata clear

// RUN FILES NEEDED FOR ESTIMATION // 
clear
run helpfunctest.ado
run estimates.ado
run classify.ado
use policy_rate.dta
oprobit y bias house gdp spread
predict p1 p2 p3
classify y, var(p1 p2 p3)