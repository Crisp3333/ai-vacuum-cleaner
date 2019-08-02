% Course: CS7075
% Student name: Dane Hylton
% Student ID: 000-74-0557
% Assignment #3 real world perceptron: 
% Due Date: 12/02/20016
% Signature: ______________
% Score: ______________

n = 6;
world = round(rand(n,n));
a = NaN(1,n);
worldChange = world;
worldChange = worldChange.';
worldChange = [worldChange; a];
a = [a NaN];
a = a.';

worldChange = [worldChange a];


% generate random -1 or 1
x1 = 0; 
x2 = 0;
x3 = 0; % unbiased input
x4 = 0;
x5 = 0; % control

sx1 = 0.49;
sx2 = 0.54;
sx3 = 0.50;
sx4 = 0.52;
sx5 = -0.5;

soutx1 = 0.5;
soutx2 = 0.1;
soutx3 = -0.3;
soutx4 = 0.8;
soutx5 = -0.1;

soutx = 0;

%round(1.49999999)
%A = ones(1) - floor(rand(1)*2)*2

% weights
w1 = (rand*2) - 1;
w2 = (rand*2) - 1;
w3 = (rand*2) - 1;
w4 = (rand*2) - 1;
w5 = (rand*2) - 1;

% Error value
delta1 = 0;
delta2 = 0;
delta3 = 0;
delta4 = 0;
delta5 = 0;

rate = 0.05; % learning rate;

% setting SS
sout = 0; % 1 when it is the correct output
pout = 0;
ss = -0.1;
tE = -1;

% Outputs
oE = 0; % output
tE = 0; % Intended output

% Network intended outputs

x1_left = 1;
x2_right = 2;
x3_up = 3;
x4_down = 4;
control = 5;
direction = 0;
% To collect weights
percept = [];
loops = 0; % number of iteration
iter = 0;
x = 1;
y = 1;
pos = [x y];
prev = [];
currentx = 0;
currenty = 0;
cleanWorld = world;
last_direction = 0;
loop = 0;
bestvalue = [];

posvec = [];
prev = [x y];
train = [];
dispWorld = [];

dirt = imread('redirt.jpg');
vacuum = imread('vacum.jpg');
posvex = [];
while(iter < 41)
    
    
    x1 = -1; 
    x2 = -1;
    x3 = -1; % unbiased input
    x4 = -1;
    x5 = 1; % control
    
    newpos = 0;
    pos1 = [];
    pos2 = [];
    pos3 = [];
    pos4 = [];
    
    posx1 = [];
    posx2 = [];
    posx3 = [];
    posx4 = [];
    
    %posx1
    direction = 0;
    %x = currentx;
    %y = currenty;
    
    loops = 0; % Reset loops
    
    %assigning input if dirt is a next tile
    if ((x - 1) ~= 0) % Move left
        xup = 1;
        if (cleanWorld(x-1, y) == 1)
            currentx = x - 1; 
            pos1 = [currentx, y];
            x1 = 1;
        else
            currentx = x - 1; 
            posx1 = [currentx, y];
        end 
    end
    if((x + 1) ~= (n + 1)) % Move right   
        if (cleanWorld(x + 1, y) == 1)
            currentx = x + 1; 
            pos2 = [currentx, y];
            x2 = 1;
        else
            currentx = x + 1; 
            posx2 = [currentx, y];
        end 
    end
    if((y - 1) ~= 0) % Move up   
        if (cleanWorld(x, y - 1) == 1)
            currenty = y - 1; 
            pos3 = [x, currenty];
            x3 = 1;
        else
            currenty = y - 1; 
            posx3 = [x, currenty];
        end
    end   
    if((y + 1) ~= (n + 1)) % Move down   
        if (cleanWorld(x, y + 1) == 1)
            currenty = y + 1; 
            pos4 = [x, currenty];
            x4 = 1;
            %direction = x4_down;
        else
            currenty = y + 1; 
            posx4 = [x, currenty];
        end 
    end
    
    % Setting 
    if (x1 == 1)
        tE = 1;
        sxout = sx1;
        direction = 1;
        pos = pos1;
    elseif(x3 == 1)
        tE = 1;
        sxout = sx3;
        direction = 3;
        pos = pos3;
    elseif(x2 == 1)
        tE = 1;
        sxout = sx2;
        direction = 2;
        pos = pos2;
    elseif(x4 == 1)
        tE = 1;
        sxout = sx4;
        direction = 4;
        pos = pos4;
    else
        sxout = sx5;
        tE = -1;
    end
    
    posvec = [posvec; pos1; pos2; pos3; pos4;];
    
    
    % Inner while loop
    % Calculations
    while(newpos ~= 1)
        
        bestvalue = [];
        currentw1 = w1;
        currentw2 = w2;
        currentw3 = w3;
        currentw4 = w4;
        currentw5 = w5;

        if (loops == 0) % condition initial weights
            s = (x1 * w1) + (x2 * w2) + (x3 * w3) + (x4 * w4) + (x5 * w5);
        else % condition calculate new weights

            % Calculate Error term 
            delta1 = rate *(tE - oE) * x1;
            delta2 = rate *(tE - oE) * x2;
            delta3 = rate *(tE - oE) * x3;
            delta4 = rate *(tE - oE) * x4;
            delta5 = rate *(tE - oE) * x5;

            % Calclate new weights
            w1 = currentw1 + delta1;
            w2 = currentw2 + delta2;
            w3 = currentw3 + delta3;
            w4 = currentw4 + delta4;
            w5 = currentw5 + delta5;

            % Calculate new for threshold function
            s = (x1 * w1) + (x2 * w2) + (x3 * w3) + (x4 * w4) + (x5 * w5);
        end
        
        % collect data
        
        % Network outputs
        % Check the output condition
        if( tE == 1) % condition with at least a 1
           if(s > sxout)
                sout = 1;
                oE = 1;
                newpos = 1;
            else
                newpos = 2;
                oE = -1;
            end
        else
            if(s <= sxout)
                sout = 1;
                oE = -1;
                
            else
                oE = 1;
            end
            newpos = 0;
        end
        
        % Vacuum has found dirt
        if (newpos == 1)
            x = pos(1,1);
            y = pos(1,2);
            cleanWorld(x, y) = 0; % clean cell
            prev = [prev; pos];
            
            last_direction = direction;
            for ii = 1:n
                for jj=1:n
                    worldChange(ii,jj) = cleanWorld(ii,jj);
                end
            end
            clf
            % for displaying the vacuum
            delete(findall(gcf,'Tag','8puzzle'))
            posit=get(gca,'position');

            [xxx , yyy] = size(cleanWorld);
            width=posit(3)/(yyy);
            height =posit(4)/(xxx);
            Cmap = [1 1 1];
            colormap(Cmap);
            pcolor(worldChange)
            axis off
            for i=1:xxx
                for j=1:yyy
                    if(cleanWorld(i,j) == 1)
                        %axes off
                        axes('posit',[posit(1)+width*(i-1),posit(2)+height*(j-1),width,height], 'Tag','8puzzle');

                        imshow(dirt)
                        %imshow('dirt.jpg');
                    end
                    if (i == x && j == y)
                        axes('posit',[posit(1)+width*(i-1),posit(2)+height*(j-1),width,height], 'Tag','8puzzle');
                        imshow(vacuum)
                    end
                end
            end
            

            set(gca, 'Ydir', 'reverse');
            pause(1)
            
            
        % Else condition here    
        elseif (newpos == 0)
            tE = 1;           
            %bestvalue = [abs(s-sx1) abs(s-sx2) abs(s-sx3) abs(s-sx4)];
            if (last_direction == 1)
                if(isempty(posx1) == 0)
                    bestvalue = [bestvalue; abs(s - w1*sx1)];
                end
                
                if(isempty(posx3) == 0)
                    bestvalue = [bestvalue; abs(s - w3*sx3)];
                end
                if(isempty(posx4) == 0)
                    bestvalue = [bestvalue; abs(s - w4*sx4)];
                end
                 
                % giving direction
                [M,I] = min(bestvalue(:)); % find smallest index
                
                if(M == abs(s - w1*sx1))
                    x1 = 1;
                    sxout = w1;
                    pos = posx1;
                    direction = 1;
                elseif(M == abs(s - w3*sx3))
                    x3 = 1;
                    sxout = sx3;
                    pos = posx3;
                    direction = 3;
                else
                    x4 = 1;
                    sxout = sx4;
                    pos = posx4;
                    direction = 4;
                end
                %pos                
                
            elseif(last_direction == 2)             
                
                if(isempty(posx2) == 0)
                    bestvalue = [bestvalue; abs(s - w2*sx2)];
                end
                if(isempty(posx3) == 0)
                    bestvalue = [bestvalue; abs(s - w3*sx3)];
                end
                
                if(isempty(posx4) == 0)
                    bestvalue = [bestvalue; abs(s - w4*sx4)];
                end
                         
                [M,I] = min(bestvalue(:)); % find smallest index
                
                if(M == abs(s - w2*sx2))
                    x2 = 1;
                    sxout = sx2;
                    pos = posx2;
                    direction = 2;
                elseif(M == abs(s - w3*sx3))
                    x3 = 1;
                    sxout = sx3;
                    pos = posx3;
                    direction = 3;
                else
                    x4 = 1;
                    sxout = sx4;
                    pos = posx4;
                    direction = 4;
                end
                
            elseif(last_direction == 3)
                
                if(isempty(posx1) == 0)
                    bestvalue = [bestvalue; abs(s - w1*sx1)];
                end
                
                if(isempty(posx2) == 0)
                    bestvalue = [bestvalue; abs(s - w2*sx2)];
                end
                
                if(isempty(posx3) == 0)
                    bestvalue = [bestvalue; abs(s - w3*sx3)];
                end               
                
                [M,I] = min(bestvalue(:)); % find smallest index
                
                if(M == abs(s - w1*sx1))
                    x1 = 1;
                    sxout = sx1;
                    pos = posx1;
                    direction = 1;
                elseif(M == abs(s - w2*sx2))
                    x2 = 1;
                    sxout = sx2;
                    pos = posx2;
                    direction = 2;
                else
                    x3 = 1;
                    sxout = sx3;
                    pos = posx3;
                    direction = 3;
                end
                               
            elseif(last_direction == 4)
                %bestvalue = [abs(s - sx1) abs(s - sx2) abs(s - sx3)];
                
                if(isempty(posx1) == 0)
                    bestvalue = [bestvalue; abs(s - w1*sx1)];
                end
                
                if(isempty(posx2) == 0)
                    bestvalue = [bestvalue; abs(s - w2*sx2)];
                end
                
                if(isempty(posx4) == 0)
                    bestvalue = [bestvalue; abs(s - w4*sx4)];
                end
                
                [M,I] = min(bestvalue(:)); % find smallest index
                
                if(M == abs(s - w1*sx1))
                    x1 = 1;
                    sxout = sx1;
                    pos = posx1;
                    direction = 1;
                elseif(M == abs(s - w2*sx2))
                    x2 = 1;
                    sxout = sx2;
                    pos = posx2;
                    direction = 2;
                else
                    x4 = 1;
                    sxout = sx4;
                    pos = posx4;
                    direction = 4;
                end
                
            % This case no dirt in any moveable positions    
            else % Last_position = 0 initial position 
                if(isempty(posx1) == 0)
                    bestvalue = [bestvalue; abs(s - w1*sx1)];
                end
                
                if(isempty(posx2) == 0)
                    bestvalue = [bestvalue; abs(s - w2*sx2)];
                end
                
                if(isempty(posx3) == 0)
                    bestvalue = [bestvalue; abs(s - w3*sx3)];
                end
                if(isempty(posx4) == 0)
                    bestvalue = [bestvalue; abs(s - w4*sx4)];
                end
                
                % Choosing the best value
                [M,I] = min(bestvalue(:)); % find smallest index
                
                if(M == abs(s - w1*sx1))
                    x1 = 1;
                    sxout = sx1;
                    pos = posx1;
                    direction = 1;
                    
                elseif(M == abs(s - w2*sx2))
                    x2 = 1;
                    sxout = sx2;
                    pos = posx2;
                    direction = 2;
                elseif(M == abs(s - w3*sx3))
                    x3 = 1;
                    sxout = sx3;
                    pos = posx3;
                    direction = 3;
                else
                    x4 = 1;
                    sxout = sx4;
                    pos = posx4;
                    direction = 4;                   
                end
                
            end % End of choosing input            
        end
        posvex = [pos last_direction];
        train = [train; x1 x2 x3 x4 x y s sxout last_direction];
        
        
        %percept = [percept; x1 x2 x3 w1 w2 w3 tE s sxout loops];
        %disp(percept)
        loops = loops + 1;
    
    end
     
    %pause(0.1)
    iter = iter + 1;
    %disp(cleanWorld)
end
%disp(percept)