// SPECIFY TO THE LOCATION OF YOUR FILES // 
cd "C:\Users\janwi\OneDrive\Documents\Classify\Classify\Old Swopit commands" 
mata: mata clear

// RUN FILES NEEDED FOR ESTIMATION // 
run helpfunctest.ado
run estimates.ado
run classify.ado

classify, var(testvar1 testvar2)