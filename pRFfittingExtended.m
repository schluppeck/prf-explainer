%% An example of pRF model fitting 
%
% this script explores using regularised linear model solutions 
% (lasso/ridge) - as provided by matlab stats toolbox.
%
% Denis Schluppeck, University of Nottingham
% 
% <https://twitter.com/schluppeck @schluppeck>
% 
% 2019-06-05

%% Background
% 
% The Dumolin & Wandell ("direct fitting") method for estimating
% population receptive fields relies on fitting a reduced, non-linear model
% to the data. There are various elaborations of the initial model (simple
% 2d gaussian, with one std parameter $[x_0, y_0, \sigma]$.
%
% A related method that doesn't make assumptions about pRF shape was
% presented in a paper from the Smirnakis group:
% 
%   http://smirnakislab.bwh.harvard.edu/wp-content/uploads/2017/01/Lee2013by_NeuroImage_2013-1.pdf
% 
% The method relies on linear approach with some regularisation to deal
% with the fact that (# of parameters) >> (# of data points). There are
% lots of different options, including rigde regression / lasso / svm
% regression, etc.
%
%

%% Load in an example stimulus. Information about the stimulus is stored in 
% this |struct|. The |im| field essentially contains a little video of the stimulus 
% aperture (white, where there was a stimulus; black where not).

load('pRFStimImage_example.mat');
stimImage = pRFStimImage;
% make sure all data are "double"
stimImage.im = double(stimImage.im);

%% visualise
% and look at it unwrapped in time (as a |montage|). Make sure to permute 
% x and y...

figure
montage(permute(stimImage.im, [2,1,3]), ...
    'BorderSize', [2,2], ...
    'BackgroundColor',[1,1,1]*0.5)


%% Example of a population receptive field (pRF)
% The following makes use of some functions in the |mgl| toolbox, so make sure 
% this is on the path. What's on the x| and |y| |axes can be a bit confusing, 
% but essentially this is a description of the area of visual space that we think 
% a particular voxel cares about ? the *population receptive field. This is* interactive!

% vectors
x = stimImage.x(:,1,1);
y = squeeze(stimImage.y(1,:,1));

% mesh versions
mx = stimImage.x;
my = stimImage.y;


%% How would a pRF respond to the stimulus? (Dumoulin & Wandell, direct method)
%
% A loop for calculating the pRF response at each point in time can be
% avoided by the use of a bit of linear algebra. Only a bit trickier to grok 
% than a loop? but quite a lot faster to calculate. 
%
% Reshape the spatial dimension of the stimulus into one big long vector (ie. 
% string all the pixels in the stimulus into one big long thing of length |nx*ny|), 
% do the same for the |pRFhat| and use a simple dot produce (linear algebra) to 
% calculate the response.
% 
% $$\mathbf{R_t} = \mathbf{\text{stimulus}_{xy,t}^{T}}\cdot\mathbf{\text{pRF}_{xy}}$$

%% direct method (forward equation) for the pRF response w/o HRF
% 
%   R(t) = reshape(stimulus, nx*ny, nt)' * pRFhat
%

%% Thinking about haemodynamics
% next step: convolve this response with an HRF shape to get an estimate of 
% an actual fMRI response that we'd expect to see. The function |calculate_prf_response() 
% |is written in such a way that you provide an HRF as an additional input argument, 
% it returns the HRF-convolved version.

% the time between image acquisitions:
TR = mean(diff(stimImage.t));
t_hrf = 0:TR:30;
hrf = fmribHRF(t_hrf);

% plot(t_hrf, hrf, 'b', 'linewidth',2)
% xlabel('Time (s)');
% ylabel('HRF impulse response')

%% 
% With this shape of the haemodynamics (assumed for now, but could also be 
% made part of the fitting), we can now calculate the expected fMRI response:


%% 
% A closure (anonymous function) that calculates response based on a simple 
% set of parameters... this is a function of only the pRF shape:
%

calculate_with_p = @(p) p(4).*calculate_prf_response(stimImage, ...
    mgauss(p, mx, my), ...
    'reshape', hrf);

%% What about actual data?
% To find out what the best parameter choices| theta =[x0, y0, sigma_x, sigma_y] 
% |are, we need to optimize (minimize the error between a model prediction for 
% a particular choice of |theta| and the data.
% 
% Let's load in a timseries (from scan #8, voxel |[10, 31, 18]| in the dataset 
% from JG's excellent tutorial on his implementation: <http://gru.stanford.edu/doku.php/mrTools/tutorialsprf) 
% http://gru.stanford.edu/doku.php/mrTools/tutorialsprf)>
%
load('tSeries_10_31_18.mat', 'tSeries')


%% Putting this together:
% Now, let's see how changes in the pRF parameter lead to differences in the 
% functional imaging response ? more directly without the intermediate visualisations. 

%% Actual fitting
% Dragging sliders by hand makes for a nice demonstration of how parameter choices 
% chang the model output. The real aim of the analysis, however, is to find the 
% optimal parameter choices pHat, that minimise some loss between the model and 
% data. Least squares (the L2 norm) is a good and commonly used choice.
% 
% There are many ways to do the parameter searching, optimisation. But an 
% easy to understand & teach way is to use |lsqnonlin().| Adapting the example 
% in the help for that function:

%%
ydata = tSeries;             % example ydata     
p0 = [-3.0, -3.0, 1.0, 10000]; % starting vals
fitFunc = @(p) calculate_with_p(p) - ydata;
lb = [-inf, -inf, 0, 1];
ub = +inf(1,4);

pHat = lsqnonlin(fitFunc, p0, lb, ub);
pHat'

%% 
% Now feed those parameters back into our function to see what the best 
% fitting model shape is:
%
tSeriesHat = calculate_with_p(pHat);

figure
subplot(3,1,[1 2])
p_ = plot(stimImage.t, tSeries, 'k-', stimImage.t, tSeriesHat, 'r-');
set(p_, 'LineWidth',2)
% the axis limits
axLimits = cell2mat(get(gca, {'xlim', 'ylim'}));
title(sprintf('Data and best fit w/ p=[%.2f, %.2f, %2.f, %.0f]',pHat))
 
% The residuals are simply the difference between model and data:
subplot(3,1,3)
plot(stimImage.t, tSeries - tSeriesHat, 'b', 'LineWidth',2);
% use the same axis limits
axis(axLimits)
title('Residuals (data-model)')

%% Estimated pRF

% calculate
pRFhat = mgauss(pHat, mx, my);

% display
figure
imagesc(x,y,permute(pRFhat, [2 1]))
hold('on')
line([0,0], get(gca,'ylim'),'color','w', 'linewidth',2);
line(get(gca,'xlim'),[0,0],'color','w', 'linewidth',2);
axis('image')
colormap(parula())
xlabel('Visual space (x)');
ylabel('Visual space (y)');

%% Lee et al method (lasso / ridge)

% given |hrf| and reshaped stimImage.im (# of pixels in stimspace by #
% timepoints)
%
% the setup that several of these methods share is that they pose the
% problem as a ill-constrained linear model
%
% eqs 2/3 in paper
% http://smirnakislab.bwh.harvard.edu/wp-content/uploads/2017/01/Lee2013by_NeuroImage_2013-1.pdf
%


%% the setup of the analysis, then, is

% Response = K matrix * weights
%
% R = K * b
%
% R is measured, K is defined above, desired: best **b**

K = estimate_prf_linear_transform(pRFStimImage,hrf);

% make sure random seeds is is set; make reproducible 
rng(42)

tic
Blasso = fitrlinear(K, tSeries, 'Regularization', 'lasso' );
Bridge =  fitrlinear(K, tSeries, 'Regularization', 'ridge' );
Brsvm =  fitrsvm(K, tSeries);
toc

%% function for reshaping vector -> image

reshape_estimate = @(x) reshape(x, size(mx));

estimatedPRFlasso = reshape_estimate(Blasso.Beta);
estimatedPRFridge = reshape_estimate(Bridge.Beta);
estimatedPRFsvm = reshape_estimate(Brsvm.Beta);



%% and now use to calculate model prediction by 
%  substituting back in

lassoOut = calculate_prf_response(stimImage,...
    estimatedPRFlasso, 'reshape', hrf);

ridgeOut = calculate_prf_response(stimImage,...
    estimatedPRFridge, 'reshape', hrf);

svmOut = calculate_prf_response(stimImage,...
    estimatedPRFsvm, 'reshape', hrf);

% https://stats.stackexchange.com/questions/48045/can-the-bias-introduced-by-lasso-change-the-sign-of-a-coefficient
% http://statweb.stanford.edu/~tibs/lasso/lasso.pdf
% check sign of correlation of OLS and LASSO

signForLasso = sign( corr(lassoOut, tSeriesHat) );
signForRidge = sign( corr(ridgeOut, tSeriesHat) );
signForSVM = sign( corr(svmOut, tSeriesHat) );

% display
figure
subplot(4,1,1)
imagesc(x,y,permute(1.0 .* pRFhat, [2 1]))
formatPlot()

subplot(4,1,2)
imagesc(x,y,permute(signForLasso .* estimatedPRFlasso, [2,1]))
formatPlot()

subplot(4,1,3)
imagesc(x,y,permute(signForRidge .* estimatedPRFridge, [2,1]))
formatPlot()


subplot(4,1,4)
imagesc(x,y,permute(signForSVM .* estimatedPRFsvm, [2,1]))
formatPlot()


%%


figure

dataZ = zscore(tSeries);

directFitZ = zscore(tSeriesHat);
lassoZ = zscore(signForLasso .* lassoOut);
ridgeZ = zscore(signForRidge .* ridgeOut);
svmZ = zscore(signForSVM .* svmOut);

subplot(2,1,1)
p_ = plot(stimImage.t, dataZ, 'k-', ...
    stimImage.t, directFitZ , 'r-');
set(p_, 'LineWidth',2)

title(sprintf('Data and best fit w/ p=[%.2f, %.2f, %2.f, %.0f]',pHat))
legend(p_, {'data', 'direct fit'})
 
% The residuals are simply the difference between model and data:
subplot(2,1,2)
p2_ = plot(stimImage.t, dataZ, 'k-', ...
    stimImage.t, lassoZ, 'm', ...
    stimImage.t, ridgeZ, 'r', ...
    stimImage.t, svmZ, 'b');
set(p2_, 'LineWidth',2);
legend(p2_, {'data',  'lasso', 'ridge', 'svm'})


%% cross-validation !
%
%
% if fitrlinear is run with 'crossvalidation' the the object that gets
% returned in a bit more complicated.

% if strcmpi(class(Blasso), 'classreg.learning.partition.RegressionPartitionedLinear')
%     % cross-validated
%     estimatedPRFlasso = reshape(Blasso.Trained{10}.Beta, size(mx));
%     
%     mse = kfoldLoss(Blasso);
%     
% else
%     %one shot
%     estimatedPRFlasso = reshape(Blasso.Beta, size(mx));
%     
% end



%% Notes and references
% Tested with |MATLAB Version: 9.4.0.813654 (R2018a)| on macOS. If you find 
% errors, typos, better ways of doing things, <https://github.com/schluppeck/prf-explainer/issues 
% submit and issue on the github repo>.
% 
% 
% 
% * Serge Dumoulin's homepage (original paper, software) ? <http://www.spinozacentre.nl/dumoulin/ 
% http://www.spinozacentre.nl/dumoulin/>
% * Brian Wandell's homepage ? <https://web.stanford.edu/group/vista/cgi-bin/wandell/ 
% https://web.stanford.edu/group/vista/cgi-bin/wandell/>
% * Justin Gardner's homepage ? <http://gru.stanford.edu/doku.php/shared/home 
% http://gru.stanford.edu/doku.php/shared/home>
% * <http://gru.stanford.edu/doku.php/mrTools/tutorialsprf http://gru.stanford.edu/doku.php/mrTools/tutorialsprf> 
% ? mrTools implementation (and tutorial for how to use)
% 

%% helper function
function formatPlot()
% formatPlot - format plot and add v/h lines

hold('on')
line([0,0], get(gca,'ylim'),'color','w', 'linewidth',2);
line(get(gca,'xlim'),[0,0],'color','w', 'linewidth',2);
axis('image')
colormap(parula())
colorbar()
xlabel('Visual space (x)');
ylabel('Visual space (y)');

end

