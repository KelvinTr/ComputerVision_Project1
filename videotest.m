clear;
I=imread('test_image.bmp');

[lx, ly]=size(I);
%figure(1), imshow(I);
load observe.dat
[On, Ot]=size(observe);
for i=1:On
    mx=observe(i,1);
    my=observe(i,2);
    for j=mx-2:mx+2
        for k=my-2:my+2
            I(k,j)=Ot;
        end
    end
end
figure(2), imshow(I)
