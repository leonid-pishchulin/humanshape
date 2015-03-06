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

function c = eval_get_color_new(cidxs)

% http://colorbrewer2.org/

color = cell(7,1);

% qualitative
color{1}{1} = [55,126,184];
color{2}{1} = [255,127,0];
color{3}{1} = [255,255,51];%[178,223,138];
color{4}{1} = [0,0,0];
color{5}{1} = [77,175,74];
color{6}{1} = [228,26,28];
% color{6}{1} = [255,255,51];
color{7}{1} = [152,78,163];
color{8}{1} = [247,129,191];
color{9}{1} = [216,179,101];

% sequential
color{1}{2} = [253,141,60];
color{1}{3} = [254,204,92];
color{1}{4} = [255,255,178];
color{1}{5} = [254,240,217];

color{2}{2} = [107,174,214];
color{2}{3} = [189,215,231];
color{2}{4} = [117,107,177];

color{3}{2} = [194,230,153];

color{4}{2} = [140,150,198];
color{4}{3} = [179,205,227];
color{4}{4} = [237,248,251];

c = color{cidxs(1)}{cidxs(2)};

end