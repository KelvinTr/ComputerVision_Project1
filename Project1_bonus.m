clear;

% Image, 2D Points, and 3D Points Loading
I1=imread('rig.bmp');
V1=imread('rig.bmp');

px = [248 239; 261 225; 270 215; 281 203; 289 194; 300 184; 307 173; 316 163; 323 155];
py = [232 234; 215 228; 199 223; 183 217; 168 211; 154 207; 141 202; 127 197];
pz = [243 205; 239 169; 234 128; 228 86; 298 30; 317 15];
pxz = [258 210; 278 187; 296 167; 314 148; 254 175; 275 154; 294 139; 312 114; 249 136; 271 115; 292 96; 331 79; 244 96; 267 76; 289 58; 308 41];
pyz = [228 217; 197 207; 165 197; 138 187; 223 182; 190 172; 159 164; 129 155; 217 145; 182 136; 150 127; 120 119; 210 103; 176 95; 142 89; 112 82];
p = [px; py; pz; pxz; pyz];

Px =  [0 0 0 1; 1 0 0 1; 2 0 0 1; 3 0 0 1; 4 0 0 1; 5 0 0 1; 6 0 0 1; 7 0 0 1; 8 0 0 1];
Py =           [0 1 0 1; 0 2 0 1; 0 3 0 1; 0 4 0 1; 0 5 0 1; 0 6 0 1; 0 7 0 1; 0 8 0 1];
Pz =           [0 0 2 1; 0 0 4 1; 0 0 6 1; 0 0 8 1; 7 0 8 1; 8 0 8 1];
Pxz =          [1 0 1 1; 3 0 1 1; 5 0 1 1; 7 0 1 1; 1 0 3 1; 3 0 3 1; 5 0 3 1; 7 0 3 1; 1 0 5 1; 3 0 5 1; 5 0 5 1; 7 0 5 1; 1 0 7 1; 3 0 7 1; 5 0 7 1; 7 0 7 1];
Pyz =          [0 1 1 1; 0 3 1 1; 0 5 1 1; 0 7 1 1; 0 1 3 1; 0 3 3 1; 0 5 3 1; 0 7 3 1; 0 1 5 1; 0 3 5 1; 0 5 5 1; 0 7 5 1; 0 1 7 1; 0 3 7 1; 0 5 7 1; 0 7 7 1];
P = [Px; Py; Pz; Pxz; Pyz];

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
[V D] = eig(Q'*Q);
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
rho = 1/norm(a3, 1);
num = -(dot(cross(a1,a3), cross(a2,a3)));
den = norm(cross(a1,a3),1)*norm(cross(a2,a3),1);

u0 = rho*rho*dot(a1,a3)                             % u0
v0 = rho*rho*dot(a2,a3)                             % v0
theta = acos(num/den)                               % theta
alpha = rho*rho*norm(cross(a1,a3),1)*sin(theta)     % alpha
beta = rho*rho*norm(cross(a2,a3),1)*sin(theta)      % beta

% Extrinsic Parameters - Rotation Matrix, and Shift Matrix
r1 = cross(a2,a3)/norm(cross(a2,a3));
r3 = rho*a3;
r2 = cross(r3, r1);

R = [r1;r2;r3]                                      % Rotation Matrix

K = [alpha, -1*alpha*cot(theta), u0; 0, beta/sin(theta), v0; 0, 0, 1]';
t = rho*K\b                                         % Shift Matrix

%% Calculate Points with Estimated Projection Matrix and Parameters
for z=1:9      %{x,10,z|x,z=0,…,10}
    for y=1:9
        point = [(M*[0, y-1, z-1, 1]')];
        point = round((point/point(3,:))');
        if y==1 && z==1
            p2 = point(:,1:2);
        else
            p2 = [p2; point(:,1:2)];
        end
    end
end

for x=1:9      %{10,y,z|y,z=0,…,10}
    for z=1:9
        point = [(M*[x-1, 0, z-1, 1]')];
        point = round((point/point(3,:))');
        if z==1 && x==1
            p3 = point(:,1:2);
        else
            p3 = [p3; point(:,1:2)];
        end
    end
end

%% Plot the estimated 2D-points onto image
for i=1:size(p2,1)          %{x,10,z|x,z=0,…,10}
    mx=p2(i,1);
    my=p2(i,2);
    for j=mx-2:mx+2
        for k=my-2:my+2
            I1(k,j)=1;
        end
    end
end
for i=1:size(p2,1)          %{10,y,z|y,z=0,…,10}
    mx=p3(i,1);
    my=p3(i,2);
    for j=mx-2:mx+2
        for k=my-2:my+2
            I1(k,j)=1;
        end
    end
end
%figure(1), imshow(I1)

%% Video - Moving 3D Object
vidObj = VideoWriter("3D_Movement_Bonus.avi");
open(vidObj);
cube_pos = [0,-1,0; 0,0,0; 1,-1,0; 1,0,0; 0,-1,1; 0,0,1; 1,-1,1; 1,0,1];
for k=1:740
    % Estimate 2D coordinate for current 3D position
    for c=1:size(cube_pos,1)
        pos_increment = [cube_pos(c,:), 1];
        new_pos = [(M*pos_increment')];
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
    if k < 141      % (0,0,0) -> (5,0,0)
        cube_pos = cube_pos + [1/20,0,0; 1/20,0,0; 1/20,0,0; 1/20,0,0; 1/20,0,0; 1/20,0,0; 1/20,0,0; 1/20,0,0];
    elseif k < 201  % (10,5,0) -> (10,5,5)
        cube_pos = cube_pos + [0,0,1/20; 0,0,1/20; 0,0,1/20; 0,0,1/20; 0,0,1/20; 0,0,1/20; 0,0,1/20; 0,0,1/20];
    elseif k < 361  % (5,5,0) -> (10,5,0)
        cube_pos = cube_pos + [-1/20,0,0; -1/20,0,0; -1/20,0,0; -1/20,0,0; -1/20,0,0; -1/20,0,0; -1/20,0,0; -1/20,0,0];
    elseif k < 521  % (10,5,5) -> (10,10,5)
        cube_pos = cube_pos + [0,1/20,0; 0,1/20,0; 0,1/20,0; 0,1/20,0; 0,1/20,0; 0,1/20,0; 0,1/20,0; 0,1/20,0];
    elseif k < 601  % (10,5,0) -> (10,5,5)
        cube_pos = cube_pos + [0,0,1/20; 0,0,1/20; 0,0,1/20; 0,0,1/20; 0,0,1/20; 0,0,1/20; 0,0,1/20; 0,0,1/20];
    else            % (10,10,5) -> (10,10,10)
        cube_pos = cube_pos + [0,-1/20,0; 0,-1/20,0; 0,-1/20,0; 0,-1/20,0; 0,-1/20,0; 0,-1/20,0; 0,-1/20,0; 0,-1/20,0];
    end
end
close(vidObj);
