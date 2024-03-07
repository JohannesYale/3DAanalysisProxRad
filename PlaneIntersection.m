% Copyright (c) 2024 Johannes Sieberer
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
function [intersectionPoints] = PlaneIntersection(mesh, intHeight)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
vertices = mesh.vertices;
faces = mesh.faces;
intersectionPoints =[];
for m = 1: height(faces)
vert1 =  vertices(faces(m,1),:);
vert2 =  vertices(faces(m,2),:);
vert3 =  vertices(faces(m,3),:);

    if(abs(sign(vert1(3)-intHeight)+sign(vert2(3)-intHeight)+sign(vert3(3)-intHeight))==1)
        intersectionDoPoints = [];
           if(sign(vert1(3)-intHeight)+sign(vert2(3)-intHeight)==0)
            intersectionVert = vert2 - vert1;
            intersectionVertScale =(vert1(3)- intHeight)/intersectionVert(3);
            intersectionDoPoints = [intersectionDoPoints; vert1 + intersectionVert * intersectionVertScale];
           end
           if(sign(vert1(3)-intHeight)+sign(vert3(3)-intHeight)==0)
             intersectionVert = vert3 - vert1;
             intersectionVertScale= (vert1(3)- intHeight)/intersectionVert(3);
             intersectionDoPoints = [intersectionDoPoints; vert1 + intersectionVert * intersectionVertScale];
           end
           if(sign(vert2(3)-intHeight)+sign(vert3(3)-intHeight)==0)
              intersectionVert = vert3 - vert2;
              intersectionVertScale =(vert2(3)- intHeight)/intersectionVert(3);
              intersectionDoPoints = [intersectionDoPoints; vert2 + intersectionVert * intersectionVertScale];
           end
        intersecMidPoint = (intersectionDoPoints(1,:) + intersectionDoPoints(2,:))/2;
        intersectionPoints= [intersectionPoints;intersecMidPoint];
    end
end

