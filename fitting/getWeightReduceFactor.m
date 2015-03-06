%{
    This file is part of the evaluation of the 3D human shape model as described in the paper:

    Leonid Pishchulin, Stefanie Wuhrer, Thomas Helten, Christian Theobalt and Bernt Schiele
    Building Statistical Shape Spaces for 3D Human Modeling
    ArXiv, March 2015

    Please cite the paper if you are using this code in your work.
    
    Author: Leonid Pishchulin.

    The code may be used free of charge for non-commercial and
    educational purposes, the only requirement is that this text is
    preserved within the derivative work. For any other purpose you
    must contact the authors for permission. This code may not be
    redistributed without permission from the authors.
%}

function [wrSmoothBody,wrSmoothHand,wrSmoothHead,wrLandmarks] = getWeightReduceFactor(scheduleIdx)

switch scheduleIdx
    
    case 0
        wrSmoothBody = 1;
        wrSmoothHand = 1;
        wrSmoothHead = 1;
        wrLandmarks = 1;
        
    case 1
        wrSmoothBody = 0.5;
        wrSmoothHand = 0.5;
        wrSmoothHead = 1;
        wrLandmarks = 0.5;
        
    case 2
        wrSmoothBody = 0.5;
        wrSmoothHand = 0.5;
        wrSmoothHead = 0.5;
        wrLandmarks = 0.5;
        
   case 3
        wrSmoothBody = 0.25;
        wrSmoothHand = 0.25;
        wrSmoothHead = 1;
        wrLandmarks = 0.25;
   
   case 4
        wrSmoothBody = 0.25;
        wrSmoothHand = 0.25;
        wrSmoothHead = 1;
        wrLandmarks = 1;     
        
end

end