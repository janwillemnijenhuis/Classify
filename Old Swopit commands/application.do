// SPECIFY TO THE LOCATION OF YOUR FILES // 
cd "/Users/jhuismans/Desktop/Swopit/Classify/Old Swopit commands"
mata: mata clear

// RUN FILES NEEDED FOR ESTIMATION // 
clear
run helpfunctest.ado
run estimates.ado
run classify.ado
use policy_rate.dta
gen y_gen = y + 3
oprobit y_gen bias house gdp spread
predict p1 p2 p3
classify y, var(p1 p2 p3)
classify y_gen, var(p1 p2 p3)
