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