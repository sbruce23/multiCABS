# multiCABS: Conduct conditional adaptive bayesian spectral analysis of replicated multivariate time series

Implements the procedure introduced in the paper
Zeda Li, Scott A. Bruce, Clinton J. Wutzke, and Yang Long "Conditional Adaptive Bayesian Spectral Analysis of Replicated Multivariate Time Series"
*Note that this folder contains code that implements NUT sampler, the other folders contains codes for HMC sampler.

Use DEMO.m as an illustration.

The following inputs are necessary to run the proposed method:

‘zt’ is a cell contains N replicated P dimensional time series

'uu' is a N-dimensional vector contains the covariate values 

‘varargin’ is a set of model parameters, which can be set by the program OptsMultiCABS.m  

MIT License

Copyright (c) 2020 

Permission is hereby granted, free of charge, to any person obtaining a 
copy of this software and associated documentation files (the "Software"), 
to deal in the Software without restriction, including without limitation 
the rights to use, copy, modify, merge, publish, distribute, sublicense, 
and/or sell copies of the Software, and to permit persons to whom the 
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in 
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS 
IN THE SOFTWARE.