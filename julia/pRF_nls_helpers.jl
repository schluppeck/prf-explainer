"""
    calculate_prf_response - dot product of stimImage with
    proposed pRF over many time frames 
    
    ported from matlab code with a julian twist...

    ds 2020-10-13
    
"""
function calculate_prf_response(stimImUnwrapped, pRFhat, hrf=[])
    
# idea: consider image as a 1d array!
# "vectorise" stimImage and also the pRF
    # but avoid allocation by not creating variable...?!
    # pRFhatUnwrapped = vec(pRFhat);
        
    # calculate response with one matrix product
    r = stimImUnwrapped' * vec(pRFhat);
       
    # if the user asks for the response w/ hrf included, calculate it
    # this will save time, if not!
    if isempty(hrf)
        return r
    else 
        # r_hrf = conv(r, hrf)
        r_hrf = ImageFiltering.imfilter(r, reflect(hrf)) # CONV, not CORR!
        r_hrf .-= mean(r_hrf)
        return r_hrf[1:size(stimImUnwrapped, 2)]
    end   
end
 
"""
mgauss - 2d gaussian (x,y, with params p x0,y0, sdxy)
   inputs: p (parameter vector, [x0,y0, sdxy] )
           x,y (meshgrid format)
           
          NB! params first, to make optimisation easier

ds 2019-03-20 (break out of mrTools/mgl version) 
ds 2020-10-13 julia version
"""
function mgauss(p, x, y)
    
# unpack parameters for readability
    x0 = p[1];
    y0 = p[2];
    sdxy = p[3];
    
# compute 2d gaussian
    m = exp.(-(((x .- x0).^2) ./ (2 * (sdxy.^2)) .+
           ((y .- y0).^2) ./ (2 * (sdxy.^2))));

# clamp small values to 0
    m[m .< 0.01] .= 0;

# return m
    m
end



"""
fmribHRF - hrf function used in by FMRIB / Oxford crew

Simple port of the Matlab code used elsewhere

    inputs: t, p [vector of params], scaled [default=1]
    outputs: result
    purpose: calculate a haemodynamic response function according
    to the equations given in the Jezzard et al book (page 252).
 
               H(t) = ( t / d1 ).^a1 .* exp( -(t-d1)./b1 ) 
                      .- c*( t / d2 ).^a2 .* exp( -(t-d2)./b2 ) 
 
   where di = ai*bi (default params a1, a2, b1, b2, c [6 12 0.9 0.9 0.35] )
"""
function fmribHRF(t=0.0:1.5:30.0, p=[6.0; 12.0; 0.9; 0.9; 0.35], scaled=true)

# if no parameters are passed use defaults from worsley chapter
# and a TR of 1.5s - which often used by Schluppeck et al.
#
# FMRIB book (p252)
# default params: a1, a2, b1, b2, c

# calculate the HRF
    t = collect(t);
    result = ( t ./ ( p[1] * p[3] )  ).^p[1] .* exp.(-( t .- p[1] * p[3] ) ./ p[3]) .- 
    p[5] * ( t ./ ( p[2] * p[4] )  ).^p[2] .* exp.(-( t .- p[2] * p[4] ) ./ p[4]);

    if scaled
  # normalize to SUM=1
        result = result ./ sum(result);
    end
end
