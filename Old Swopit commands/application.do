// SPECIFY TO THE LOCATION OF YOUR FILES // 
cd "/Users/jhuismans/Desktop/Swopit/Classify/Old Swopit commands"
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
