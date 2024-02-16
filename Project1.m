clear;

% Image, 2D Points, and 3D Points Loading
I1=imread('test_image.bmp');
V1=imread('test_image.bmp');
p = load("observe.dat");
P = load("model.dat");
P = [P, ones(27,1)];

%% Projection Matrix Estimation
zero_T = [0 0 0 0];
for i=1:size(p,1)
    A = [P(i,:) zero_T -1*p(i,1).*P(i,:)];
    B = [zero_T P(i,:) -1*p(i,2).*P(i,:)];
    if i == 1
        Q = [A; B];
    else
        Q = [Q; A; B];
    end
end

% Eigenvalue Problem - Evaluate which column the best error value
[V, D] = eig(Q'*Q);
min_index = 1;
for i = 1:size(D,1)
    if D(min_index, min_index) > D(i,i)
        min_index = i;
    end
end

% Projection Matrix Construction
m=V(:,min_index);
M = reshape(m,4,3)';

%% Variables to Calculate Intrinsic and Extrinsic Parameters
A = M(1:3, 1:3);
b = M(1:3, 4);
a1 = A(1,:);
a2 = A(2,:);
a3 = A(3,:);

% Intrinsic Parameters - u0, v0, theta, alpha, beta
rho = 1/norm(a3);
num = -(dot(cross(a1,a3), cross(a2,a3)));
den = norm(cross(a1,a3))*norm(cross(a2,a3));

u0 = rho*rho*dot(a1,a3)                             % u0
v0 = rho*rho*dot(a2,a3)                             % v0
theta = acos(num/den)                               % theta
alpha = rho*rho*norm(cross(a1,a3))*sin(theta)     % alpha
beta = rho*rho*norm(cross(a2,a3))*sin(theta)      % beta

% Extrinsic Parameters - Rotation Matrix, and Shift Matrix
r1 = cross(a2,a3)/norm(cross(a2,a3));
r3 = rho*a3;
r2 = cross(r3, r1);

R = [r1; r2; r3]                                      % Rotation Matrix

K = [-alpha, -1*alpha*cot(theta), u0; 0, beta/sin(theta), v0; 0, 0, 1];
t = rho*inv(K)*b                                      % Shift Matrix

M_estimated = [K*R K*t]/rho

%% Calculate Points with Estimated Projection Matrix and Parameters
for x=1:11      %{x,y,0|x,y=0,…,10}
    for y=1:11
        point = [(M_estimated*[x-1, y-1, 0, 1]')];
        point = round((point/point(3,:))');
        if y==1 && x==1
            p1 = point(:,1:2);
        else
            p1 = [p1; point(:,1:2)];
        end
    end
end

for z=1:11      %{x,10,z|x,z=0,…,10}
    for y=1:11
        point = (M_estimated*[10, y-1, z-1, 1]');
        point = round((point/point(3,:))');
        if y==1 && z==1
            p2 = point(:,1:2);
        else
            p2 = [p2; point(:,1:2)];
        end
    end
end

for x=1:11      %{10,y,z|y,z=0,…,10}
    for z=1:11
        point = (M_estimated*[x-1, 10, z-1, 1]');
        point = round((point/point(3,:))');
        if z==1 && x==1
            p3 = point(:,1:2);
        else
            p3 = [p3; point(:,1:2)];
        end
    end
end

%% Plot the estimated 2D-points onto image
for i=1:size(p1,1)          %{x,y,0|x,y=0,…,10}
    mx=p1(i,1);
    my=p1(i,2);
    for j=mx-2:mx+2
        for k=my-2:my+2
            I1(k,j)=1;
        end
    end
end
for i=1:size(p1,1)          %{x,10,z|x,z=0,…,10}
    mx=p2(i,1);
    my=p2(i,2);
    for j=mx-2:mx+2
        for k=my-2:my+2
            I1(k,j)=1;
        end
    end
end
for i=1:size(p1,1)          %{10,y,z|y,z=0,…,10}
    mx=p3(i,1);
    my=p3(i,2);
    for j=mx-2:mx+2
        for k=my-2:my+2
            I1(k,j)=1;
        end
    end
end
figure(1), imshow(I1)

%% Video - Moving 3D Object
vidObj = VideoWriter("3D_Movement.avi");
open(vidObj);
cube_pos = [0,0,0; 0,1,0; 1,0,0; 1,1,0; 0,0,1; 0,1,1; 1,0,1; 1,1,1];
for k=1:538
    % Estimate 2D coordinate for current 3D position
    for c=1:size(cube_pos,1)
        pos_increment = [cube_pos(c,:), 1];
        new_pos = (M_estimated*pos_increment');
        new_pos = round((new_pos/new_pos(3,:))');
        if c==1
            new_pos_2d = [new_pos(:,1:2)];
        else
            new_pos_2d = [new_pos_2d; new_pos(:,1:2)];
        end
    end
    
    % Cube Construction
    cube_shape = [new_pos_2d(1,:), new_pos_2d(2,:);new_pos_2d(2,:), new_pos_2d(4,:)];
    cube_shape = [cube_shape; [new_pos_2d(4,:), new_pos_2d(3,:)]];
    cube_shape = [cube_shape; [new_pos_2d(3,:), new_pos_2d(1,:)]];
    cube_shape = [cube_shape; [new_pos_2d(1,:), new_pos_2d(5,:)]];
    cube_shape = [cube_shape; [new_pos_2d(2,:), new_pos_2d(6,:)]];
    cube_shape = [cube_shape; [new_pos_2d(3,:), new_pos_2d(7,:)]];
    cube_shape = [cube_shape; [new_pos_2d(4,:), new_pos_2d(8,:)]];
    cube_shape = [cube_shape; [new_pos_2d(5,:), new_pos_2d(6,:)]];
    cube_shape = [cube_shape; [new_pos_2d(6,:), new_pos_2d(8,:)]];
    cube_shape = [cube_shape; [new_pos_2d(8,:), new_pos_2d(7,:)]];
    cube_shape = [cube_shape; [new_pos_2d(7,:), new_pos_2d(5,:)]];

    currFrame = insertShape(V1, "line", cube_shape,"LineWidth",2,"ShapeColor","cyan");
    figure(2), imshow(currFrame)
    writeVideo(vidObj, currFrame);

    % Position Increment
    if k < 101      % (0,0,0) -> (5,0,0)
        cube_pos = cube_pos + [1/20,0,0; 1/20,0,0; 1/20,0,0; 1/20,0,0; 1/20,0,0; 1/20,0,0; 1/20,0,0; 1/20,0,0];
    elseif k < 201  % (5,0,0) -> (5,5,0)
        cube_pos = cube_pos + [0,1/20,0; 0,1/20,0; 0,1/20,0; 0,1/20,0; 0,1/20,0; 0,1/20,0; 0,1/20,0; 0,1/20,0];
    elseif k < 281  % (5,5,0) -> (10,5,0)
        cube_pos = cube_pos + [1/20,0,0; 1/20,0,0; 1/20,0,0; 1/20,0,0; 1/20,0,0; 1/20,0,0; 1/20,0,0; 1/20,0,0];
    elseif k < 380  % (10,5,0) -> (10,5,5)
        cube_pos = cube_pos + [0,0,1/20; 0,0,1/20; 0,0,1/20; 0,0,1/20; 0,0,1/20; 0,0,1/20; 0,0,1/20; 0,0,1/20];
    elseif k < 460  % (10,5,5) -> (10,10,5)
        cube_pos = cube_pos + [0,1/20,0; 0,1/20,0; 0,1/20,0; 0,1/20,0; 0,1/20,0; 0,1/20,0; 0,1/20,0; 0,1/20,0];
    else            % (10,10,5) -> (10,10,10)
        cube_pos = cube_pos + [0,0,1/20; 0,0,1/20; 0,0,1/20; 0,0,1/20; 0,0,1/20; 0,0,1/20; 0,0,1/20; 0,0,1/20];
    end
end
close(vidObj);
