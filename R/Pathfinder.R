# Pathfinder

#' Generate Perturbations
#' 
#' Find curves similar to a set of example curves.  This function takes as input a set of example curves, and uses them to infer a probability distribution of curves.  \code{numPert} curves are sampled from this probability distribution.
#' @param training This is a numerical matrix of training data, where the rows represent different samples, columns represent different time points (or points on a single spatial axis), and the values correspond to measurements
#' @param tp A numerical vector of time points (or spatial coordinates along a single axis)
#' @param iterations a positive integer, representing the maximum number of iterations employed during time warping (see time_warping in fdasrvf library)
#' @param spline a positive integer, representing the degree of the B-spline interpolation when calculating values at the new, evenly spaced knot positions
#' @param knots a positive integer-- for time warping to work optimally, the points must be evenly sampled.  This determines how many points do we evenly sample before conducting time warping
#' @param numPert a positive integer, representing the number of sampled curves to output.
#' @return An fdawarp object (see fdasrvf library)
#' @examples 
#' mat=CanadianWeather$monthlyTemp
#' \donttest{generated=generatePerturbations(mat, c(1:length(mat[,1])))}

generatePerturbations<- function(training, tp, iterations=20, spline=3, knots=100, numPert=20){
    #interpolate points using splines
    newX=seq(min(tp), max(tp), (max(tp)-min(tp))/knots)
    training_interpolated=apply(training, 2, function(i){deBoorWrapper(newX, tp, i, spline)})
    #do time warping
    tw= time_warping(training_interpolated, newX, MaxItr = iterations, showplot=FALSE)
    #sample from the distribution of functions
    gm=gauss_model(tw, n=numPert)
    gm
}

# Helper function: calculates the optimisation matrix value 
optimisationMatrixValue <- function(i, j, edge, max_value, max_link, tp, y, w, multipleGenes=F, iter=20, knots=100, numPerts=1000){
if(!multipleGenes){
    #if single gene, use this method
    optimisationMatrixValueSingleGene(i, j, edge, max_value, max_link, tp, y, w)
}else{
    #take sum over the genes
    sum(sapply(c(1:dim(y)[1]), function(index){
            y1=y[index,] #check orientation
            w1_id=w[[index]]
          
                optimisationMatrixValueSingleGene(i, j, edge, max_value, max_link, tp, y1, w1_id)
                
                }))

}
    }

#Helper function: gets scoreF1
optimisationMatrixValueSingleGene <- function(i, j, edge, max_value, max_link, tp, y, w){
    #case 1: start node
    if(i==0){
        #case 1.5: start node AND end node:
        if(j>length(tp)){0}else{

        scoreF1(c(tp[1], tp), c(y[j], y), rbind(w[j,],w),  tp[1], tp[j], c(1, j+1))

        }
    }else{
        if(edge==1){
            Inf
        }else{
        index=backtrace(max_link, i, edge) 

        #case 2: end node
        if(j>length(tp)){
            if(i==length(tp)){max_value[i, edge-1]}else{
                
                scoreF1(c(tp, tp[j-1]), c(y, y[i]), rbind(w, w[i,]),  tp[i], tp[j-1], c(index, i, j))+max_value[i, edge-1]
                
            #case 3: just intermediate nodes
             }}else{
            scoreF1(tp, y, w,  tp[i], tp[j], c(index, i, j))+max_value[i, edge-1]
            }
        }
    }
   }



#' Find best subset of points for follow-up experiments, using F1 metric
#' 
#' findPathF1 finds the best subset of points to sample from a time course (or spatial axis, along a single axis), based on a set of example curves. Specifically, it finds subsets of points that estimate the shape of the curve effectively. 
#'
#' @param tp A numerical vector of time points (or spatial coordinates along a single axis)
#' @param training this is a numerical matrix of training data, where the rows represent different samples, columns represent different time points (or points on a single spatial axis), and the values correspond to measurements. (If \code{mult==TRUE}, then this is instead a list of training matrices)
#' @param numSubSamples integer that represents the number of time points that will be subsampled
#' @param spline A positive integer representing the spline used to interpolate between knots when generating perturbations.  Note that this does NOT designate the spline used when calculating the L2-error.
#' @param resampleTraining A boolean designating whether the exact training data should be used (False) or whether a probability distribution of curves should be generated and training curves resampled (True).
#' @param iter A positive integer, representing the maximum number of iterations employed during time warping (see time_warping in fdasrvf library)
#' @param knots A positive integer-- for time warping to work optimally, the points must be evenly sampled.  This determines how many points do we evenly sample before conducting time warping
#' @param numPerts a positive integer, representing the number of sampled curves to output.
#' @param fast is a boolean, which determines whether the algorithm runs in fast mode where the sum of the perturbations is calculated prior to integration.
#' @param mult is a boolean.  If mult is true, then training will be a list of training matrices.  This will be the case if there are multiple genes to consider at the same time.  Training sets will be normalised by the size of the L2-error. 
#' @param weights is a vector of numbers that is the same length as the number of training curves. This describes the relative importance of these curves.
#' @return An integer vector of the indices of the time points selected to be subsampled.  The actual time points can be found by \code{tp[output]}.  The length of this vector should be \code{numSubSamples}.
#' @examples  
#' #load data:
#' #matrix with 12 rows, representing months (time)
#' #and 35 columns, representing cities (experiments)
#' mat=CanadianWeather$monthlyTemp 
#' #find a set of points that help predict the shape of the curve:
#' \donttest{a=findPathF1(c(1:12), mat, 5, numPerts=3) #make numPerts>=20 for real data }
#' \donttest{print(a) #indices of months to select for follow-up experiments}
#' \donttest{print(rownames(CanadianWeather$monthlyTemp)[a]) #month names selected}
#' 
#' 
#' 
findPathF1 <- function(tp, training, numSubSamples, spline=1, resampleTraining=T, iter=20, knots=100, numPerts=1000, fast=T, mult=F, weights=c()){
    findPathF2(tp, rep(0, length(tp)), training, numSubSamples, spline=spline, resampleTraining=resampleTraining, iter=iter, knots=knots, numPerts=numPerts, fast=fast, mult=mult)
}

#' Find best subset of points for follow-up experiments, using F3 metric
#' 
#' findPathF3 finds the best subset of points to sample from a time course (or spatial axis, along a single axis), based on a set of example curves. Specifically, it finds subsets of points that estimate the shape of the curve, normalised by the variance. 
#'
#' @param tp A numerical vector of time points (or spatial coordinates along a single axis)
#' @param training1 this is a numerical matrix of training data of experimental condition 1, where the rows represent different samples, columns represent different time points (or points on a single spatial axis), and the values correspond to measurements.  
#' @param training2 this is a numerical matrix of training data of experimental condition 2, where the rows represent different samples, columns represent different time points (or points on a single spatial axis), and the values correspond to measurements.  
#' @param numSubSamples integer that represents the number of time points that will be subsampled
#' @param spline A positive integer representing the spline used to interpolate between knots when generating perturbations.  Note that this does NOT designate the spline used when calculating the L2-error.
#' @param resampleTraining A boolean designating whether the exact training data should be used (False) or whether a probability distribution of curves should be generated and training curves resampled (True).
#' @param iter A positive integer, representing the maximum number of iterations employed during time warping (see time_warping in fdasrvf library)
#' @param knots A positive integer-- for time warping to work optimally, the points must be evenly sampled.  This determines how many points do we evenly sample before conducting time warping
#' @param numPerts a positive integer, representing the number of sampled curves to output.
#' @param fast is a boolean, which determines whether the algorithm runs in fast mode where the sum of the perturbations is calculated prior to integration.
#' @return An integer vector of the indices of the time points selected to be subsampled.  The actual time points can be found by \code{tp[output]}.  The length of this vector should be \code{numSubSamples}.
#' 
#' @examples  
#' 
#' #Set up data:
#' namAtlantic=CanadianWeather$region[as.character(colnames(CanadianWeather$monthlyTemp))]
#' atlanticCities=which(namAtlantic=="Atlantic")
#' matAtlantic=CanadianWeather$monthlyTemp[, names(atlanticCities)]
#' 
#' namContinental=CanadianWeather$region[as.character(colnames(CanadianWeather$monthlyTemp))]
#' continentalCities=which(namContinental=="Continental")
#' matContinental=CanadianWeather$monthlyTemp[, names(continentalCities)]
#' 
#' #find a set of points that helps capture the difference 
#' #between Atlantic and Continental cities, normalised by the variance
#' #make numPerts >=20 for real data
#' \donttest{a=findPathF3(c(1:12),  matAtlantic,  matContinental, 5, numPerts=3)} 
#' \donttest{print(a) #indices of months to select for follow-up experiments}
#' \donttest{print(rownames(CanadianWeather$monthlyTemp)[a]) #month names selected}
#' 
findPathF3 <- function(tp, training1, training2, numSubSamples, spline=1, resampleTraining=F, iter=20, knots=100, numPerts=1000, fast=T){
    #generate lots of perturbations from training set 1 and 2
    generated1=generatePerturbations(training1, tp, iterations=iter, spline=spline, knots=knots, numPert=numPerts)
    generated2=generatePerturbations(training2, tp, iterations=iter, spline=spline, knots=knots, numPert=numPerts)
    
    #calculate inv coef of var
    diff=generated1$ft-generated2$ft
    variances=apply(diff, 1, function(i){var(i)})
    vals=apply(diff, 2, function(i){
        approx(generated1$time, i/sqrt(variances), xout=tp)$y
    })
   
 
    
    #plug into findPathF2
    findPathF2(tp, rep(0, length(tp)), vals, numSubSamples, spline, resampleTraining, iter, knots, numPerts, fast=fast)
}


#' Find best subset of points for follow-up experiments, using F2 metric
#' 
#' findPathF2 finds the best subset of points to sample from a time course (or spatial axis, along a single axis), based on a set of example curves. Specifically, it compares between a control curve and a set of experimental curves. 
#'
#' @param tp A numerical vector of time points (or spatial coordinates along a single axis)
#' @param y A numerical vector of measurements (of the control).  If \code{mult==TRUE}, then this will be a matrix, where each column would be the y that corresponds with each training matrix.  
#' @param training This is a numerical matrix of training data, where the rows represent different samples, columns represent different time points (or points on a single spatial axis), and the values correspond to measurements.   (If \code{mult==TRUE}, then this is instead a list of training matrices).
#' @param numSubSamples integer that represents the number of time points that will be subsampled
#' @param spline A positive integer representing the spline used to interpolate between knots when generating perturbations.  Note that this does NOT designate the spline used when calculating the L2-error.
#' @param resampleTraining A boolean designating whether the exact training data should be used (False) or whether a probability distribution of curves should be generated and training curves resampled (True).
#' @param iter A positive integer, representing the maximum number of iterations employed during time warping (see time_warping in fdasrvf library)
#' @param knots A positive integer-- for time warping to work optimally, the points must be evenly sampled.  This determines how many points do we evenly sample before conducting time warping
#' @param numPerts a positive integer, representing the number of sampled curves to output.
#' @param fast is a boolean, which determines whether the algorithm runs in fast mode where the sum of the perturbations is calculated prior to integration.
#' @param weights is a vector of numbers that is the same length as the number of training curves. This describes the relative importance of these curves.
#' @param mult is a boolean, which will determine whether multiple genes are considered at once.
#' @return An integer vector of the indices of the time points selected to be subsampled.  The actual time points can be found by \code{tp[output]}.  The length of this vector should be \code{numSubSamples}.
#' @examples 
#' #load data:
#' # a matrix with 12 rows, representing months (time) 
#' # and 35 columns, representing cities (experiments) 
#' mat=CanadianWeather$monthlyTemp 
#' y=CanadianWeather$monthlyTemp[,"Resolute"]
#' #find a set of points that help predict the shape of the curve
#' \donttest{a=findPathF2(c(1:12), y, mat, 5, numPerts=3) #make numPerts>=20 for real data}
#' \donttest{print(a) #indices of months to select for follow-up experiments}
#' \donttest{print(rownames(CanadianWeather$monthlyTemp)[a]) #month names selected}
#' 
#' 
findPathF2 <- function(tp, y, training, numSubSamples, spline=1, resampleTraining=T, iter=20, knots=100, numPerts=1000, fast=T, mult=F, weights=c()){
    
    perts=NA
    w=NA
    
    if(mult){
        
        #then training is a list of training matrices
        training=sapply(c(1:length(training)), function(index){
            
        #for each of these, generate a set of perturbations
            perts=generatePerturbations(training[[index]], tp, iterations=iter, spline=1, knots=knots, numPert=numPerts)
            a=apply(perts$ft, 2, function(i){
                deBoorWrapper(tp, perts$time, i, spline)
            })
        #calculate the sum of the area under each of these curves
            aSize=sum(sapply(c(1:numPerts), function(i){
                L2(tp, y[,index], a[,i], min(tp), max(tp), index=c(1, length(tp)))
                
            }))/numPerts
            
        #divide by the sum
            
            aAdj=sapply(c(1:numPerts), function(i){
                (a[,i]-y[,index])/aSize
                
            })
            
            #Apply weights, if necessary
            if(length(weights)==length(training)){
                aAdj=aAdj*weights[index]
            }
            
        #add up the curves 
            apply(aAdj, 1, function(i){sum(i)})
           
        })
        
        y=rep(0, length(tp))
    }
    
    
    #Usually you would want to use the fdasrvf package to generate new pdf from the training set and sample curves from that
    if(resampleTraining){
        
    
    perts=generatePerturbations(training, tp, iterations=iter, spline=spline, knots=knots, numPert=numPerts)

    w=apply(perts$ft, 2, function(i){
        deBoorWrapper(tp, perts$time, i, spline)
    })

    
    }else{
    #Sometimes a different strategy might be used to sample example curves-- 
    #in this case, training and training2 can just be set to the new set of perturbations
        w=training
    }
    if(fast){
    
        w=as.matrix(sapply(c(1:length(tp)), function(i){
           sum((w[i,]-y[i]))
        }))
        print(dim(w))
        y=rep(0, length(y))
 
   
    }
     min_score=matrix(0, nrow=1+length(tp), ncol=numSubSamples+1)
     min_link=matrix(0, nrow=1+length(tp), ncol=numSubSamples+1)

    # #Make an N x N x E
     for(j in c(1:(1+length(tp)))){
         temp=sapply(c(0:(j-1)), function(i){
             sapply(c(1:(1+numSubSamples)), function(edge){
                 if(edge>(i+1)){Inf}else{
                 optimisationMatrixValue(i, j, edge, min_score, min_link, tp, y, w, multipleGenes=mult)
                 }})
         })

         print(temp)


         min_link[j,]=apply(temp, 1, function(k){

             which.min(k)[1]-1}
            )

         min_score[j,]=apply(temp, 1, function(k){

          min(k)})


     }

     print(min_link)
     print(min_score)
     backtrace(min_link, length(min_link[,1]), length(min_link[1,]))

     ######If you allow splines for calculating the L2-error you would need to add `step 2`
     ##STEP2: do correction for first 'spline' time points that were subset, by running the algorithm backwards

}


#internal function for cal
scoreF1 <- function(tp, y, w,  start, stop, index, numSubdivisions=500){
     sum(apply(w, 2, function(i){
         tp=c(tp[1], tp, tp[length(tp)])
         y=c(y[index[1]], y, y[index[length(index)]])
         i=c(i[index[1]], i, i[index[length(index)]])
         index=c(1, index+1, length(tp))
        integrate(F1, start, stop, tp=tp, g=y, w=i, index=index, subdivisions=numSubdivisions, rel.tol=.Machine$double.eps^0.1)$value
    }))}



#' L2-error
#'
#'Given two functions y1(t) and y2(t), this function finds the L2-distance between the following two curves:
#'a) y1(t)-y2(t) sampled at all time points (\code{tp}) 
#'b) y1(t)-y2(t) sampled at the time points indexed by \code{index} (\code{tp[index]}).
#'Note that by setting \code{y2} to \code{rep(0,length(tp))}, this function can be used to estimate the L2-error in the shape of \code{y1}. 
#'
#' @param tp A numerical vector of time points (or spatial coordinates along a single axis)
#' @param y1 A numerical vector of measurements (of the control)
#' @param y2 A numerical vector of measurements (of the experimental condition)
#' @param start A numerical value representing the start time (or spatial coordinate) of the integration
#' @param stop A numerical value representing the end time (or spatial coordinate) of the integration
#' @param index A vector of positive integers representing the indices of \code{tp} that we subsample
#' @param numSubdivisions This can be adjusted to ensure the integration doesn't take too long, especially if we aren't overly concerned with rounding errors.
#'
#' @return A numeric value-- the L2 error.
L2 <-function(tp, y1, y2, start, stop, index, numSubdivisions=2000){
   tp=c(tp[1], tp, tp[length(tp)])
   y1=c(y1[index[1]], y1, y1[index[length(index)]])
   y2=c(y2[index[1]], y2, y2[index[length(index)]])
   index=c(1, index+1, length(tp))

   integrate(F1, start, stop, tp=tp, g=y1, w=y2, index=index,
               spl=1, subdivisions=numSubdivisions,rel.tol=.Machine$double.eps^0.1 )$value
}

#helper function
 meanSqr <-function(x, tp, g, w, tp2, g2, w2, spl){

    temp=deBoorWrapper(x, tp, g, spl)-deBoorWrapper(x, tp, w, spl)-(deBoorWrapper(x, tp2, g2, spl)-deBoorWrapper(x, tp2, w2, spl))
    temp*temp
}
#helper function
F1 <-function(x, tp, g, w, index, spl=1){
   # print(paste(tp[1], tp[2]))
    if(spl==1){
        if(length(index)==2 & tp[index[1]]==tp[index[2]]){
            rep(0, length(x))
        }else{
      # print(paste('tp subset', tp[2:(length(tp)-1)]))
      #      print(paste('index', index))   
      #      print(length(index)==2)
      #      print(paste('tp in ids', tp[index[1]], tp[index[2]]))
      temp=approx(tp[2:(length(tp)-1)], g[2:(length(tp)-1)]-w[2:(length(tp)-1)], xout=x)$y-approx(tp[index], g[index]-w[index], xout=x)$y 
      #print(paste('temp', temp))
      temp*temp }
    }else{
    temp=deBoorWrapper(x, tp, g, spl)-deBoorWrapper(x, tp, w, spl)-(deBoorWrapper(x, tp[index], g[index], spl)-deBoorWrapper(x, tp[index], w[index], spl))
    temp*temp
    }
    }
#helper function
deBoorWrapper <- function(x, tp, values, spline){
    if(length(tp)==2 & tp[1]==tp[2]){
        0
    }

    #check input values
    if(spline<0){
        print('spline must be >0')
    }

    sapply(x, function(i){
        k_smaller=which(tp<i)
        k=k_smaller[length(k_smaller)]

        if(i==tp[1]){values[1]}else{
        if(k<spline){
            k_smaller=which(sort(-tp)<(-i))
            k=k_smaller[length(k_smaller)]
            temp=deBoor2(k, -i, sort(-tp), values[length(values):1], spline)
            temp
        }else{
            deBoor2(k, i, tp, sapply(values, function(bl){as.numeric(bl)}), spline)
        }
        }
    })

}

#helper function
deBoor2 <- function(k, x, t, b, p){
    indices=seq(0, p, 1)+ k - p+1
    d=b[indices]

    for (r in indices[0:(length(indices)-1)]){
        temp=indices[which(indices>r)]
        for (j in temp[length(temp):1]){

            alpha = (x - t[r]) / (t[j] - t[r])
                     oldD=d[j-min(indices)+1]

            d[j-min(indices)+1] = (1.0 - alpha) * d[r-min(indices)+1] + alpha * d[j-min(indices)+1]

        }
    }
 d[p+1]
}


#Helper function: 
#finds backtrace from the max_index matrix
backtrace<- function(max_index, index, edge, until=NA){
    if(max_index[index, edge]==0){
        c();
    }else{
        if(is.na(until)){
        c(backtrace(max_index, max_index[index, edge], edge-1),max_index[index, edge])
        }else{
            if(until==0){
                c(max_index[index, edge])
            }else{
                c(backtrace(max_index, max_index[index, edge], edge-1, until=until-1),max_index[index, edge])
            }
        }
    }
}

