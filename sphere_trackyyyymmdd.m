function sphere_trackyyyymmdd
%% So, how does this code work then...
%
% The algorithm is from;
% Non-Bias-Limited Tracking of Spherical Particles, Enabling Nanometer
% Resolution at Low Magnification, Biophys. J.2362–2371 (2012)
% MTJ van Loenhout, JWJ Kerssemakers, I De Vlaminck, and C Dekker
%
% This program tracks x,y,z of an arbitrary number of microspheres
%
% It's input is 8-bit jpegs but this could probably be altered to use
% 16-bit jpegs, or tifs etc.
%
% It presumes you have the following data;
% 1. A series of images (no_steps many of them) that have been collected fom
% the microscope at discrete objective heights (1 image per objective
% height). The objective displacement in this code is treated as being '1
% vertical pixel'. In this sense you should really keep the objective step
% size constant in each experiment. These Look Up Table (LUT) images are
% named 'zstack001.jpg, zstack002.jpg,... The algorithm is designed to
% track microspheres that are out of focus. The best way to achieve this is
% focus on a ball/tether, then move the objective up to get the rings to appear.
% These images should be in a folder called yyyy-mm-dd

% 2. Your raw experimental data. I use a sequence of jpegs here named
% 'frame1.jpg,frame2.jpg,...' contained in folders named
% 'goagain0,gagain1,...'. The frame number starts again in each foler.

% You should have a zstack of the same FOV for your LUT. Each bead is
% different, lighting changes etc. So your images should be the same size
% for both raw and LUT images.
%
% The order to do things in;
%   1. set track = 0
%   2. set make_LUT = 0
%   3. Uncomment 'find positions from LUT data' part. This uses 2D cross
%   correlation to find spheres. It's output is csvlistxLUTpos.dat and csvlistyLUTpos.dat
%   in the yyyy-mm-dd_LUT folder. Basically these give the x,y coordinates of
%   all the beads
%    - run code
%   4. Use 'addpos' function to add the position of a particular singular bead
%   to the .dat files mentioned in 1.
%    - run code
%   5. Comment the previous two out when happy with spheres found and list
%   created.
%   6. set make_LUT = 1
%    - run code
%   7. in folder yyyy-mm-dd_LUT should be a file LUTposbead_noXnew.dat for
%   each sphere
%   8. set make_LUT = 0
%   9. set track = 1
%     - run code


%% run_code
directory = '';         % leave this alone - redundant

%% the code works by extracting small sqaures from the large FOV image. i.e. a little square for each sphere
roi_size = 120;         % region of interest size = square size around ball

%%
descriptor = '';  % this is just something to add to the file names

%% number of frames and folders
no_frames = 3999;        % number of frames in the folder (i.e. the number of frames ot trakc if you were tracking all of them)
start = 1;              % starting frame number
no_steps = 50;          % no of steps in LUT
folder_start = 0;       % the first folder contents you want to track
folder_end = 1;        % the last folder contents you want to track
frame_int = 100;        % here you can skip frames, to shorten process, if you want to track all images, make it 1.
frame_rate = 58;        % fps - well this is obvious, right?

%% this is the path to your data, remember raw experiment images in yyyy-mm-dd and LUT images in yyy-mm-dd_LUT
date = '//path/2050-01-01';

descriptwo = ['test' '_int' num2str(frame_int) ''];  % name for your output xyz data

%% quadrant interpolation options
no_it = 3;                  % iterations - 3 works well
no_points_z_porabola = 5;   % for sub integer z pixel precision it fits to the 5 points around the minimum difference it found.
truncate = 0;               % how many interpolted pixels to remove form the radial profile of the bead - this is the edge of the bead
cut_off = 40;%80%120;%60;   % how many interpolated pixels to remove form the start of the radial profile - this is the centre of the bead
m = 3;                      % no of interpolations. if you start with an image of the bead of roi_size x roi_size
                            % you then end up with an image that is (2*m*roi_size) x (2*m*roi_size)

%% 1 for tracking; 0 for finding balls and making LUTs
track = 1;

%% 1 for making_LUT; 0 for not
make_LUT = 0;

%% find positions from LUT data - uncomment to use
% findpositionsofbeads(directory,[date '_LUT'],'zstack010.jpg','',roi_size,0.8,0);

%% add a particular bead - uncomment to use
% addpos('',[date '_LUT'],'zstack021.jpg','',roi_size)

%% display beads with numbers - uncomment to use
% displaybeads('',[date '_LUT'],'zstack021.jpg','',roi_size)

%% read in initial positions - this reads in the positions of all the beads
xinitguess = csvread([directory date '_LUT/' 'csvlistx' descriptor '.dat'])-(roi_size/2);
yinitguess = csvread([directory date '_LUT/' 'csvlisty' descriptor '.dat'])-(roi_size/2);

%% calculate some important things
frame_ratio = floor((no_frames)/frame_int); % actual number of frames to track due to interval not being unity
no_beads = length(xinitguess);              % no. beads

startbead = 1;                              % well, this is the starting bead number...
endbead = no_beads;                         % this is the last bead to track (together these two numbers go through a list of the
% consecutive numbers 'startbead' to 'endbead'

not_to_track = linspace(1,no_beads,no_beads);   % this is the list of beads not to track - initially this says do not track
to_delete = linspace(1,no_beads,no_beads);              % write here the beads you want to track - in this case i want them all so I can use linspace
not_to_track(to_delete) = [];               % this deletes the beads you want to track from the not_to_track array

%% This whole not_to_track nonsense seems bass ackwards but it works...

%% This is Dan's special bit - you can ignore this.
%initialise bead ignore matrices
disappeared = zeros(folder_end + 1,no_beads);

%populate disappeared matrix % !!!!!! Remember that 1 = 0 for array
%indexing
wbag = (folder_end+1)*ones(1,no_beads);
% wbag(2) = 97;

if length(wbag) == no_beads
    %do nothing
else
    size(wbag);
    no_beads;
    return
end


%% LUT calculation
if make_LUT == 1
    [~] = trackerclean_9_0(directory,date,'ten',descriptor,roi_size,1,no_steps,0,no_steps,1,xinitguess,yinitguess,descriptwo,frame_rate,...
        no_it,...
        no_points_z_porabola,...
        truncate,...
        cut_off,...
        m,...
        no_steps,...
        no_beads,...
        startbead,...
        endbead,...
        not_to_track);
elseif make_LUT == 0
    % do nothing
end

%% track raw experimental data section
if track == 0
    % do nothing
elseif track == 1
    for i = folder_start:folder_end
        j = 0 + i;
        %         disp(j);
        folder_name = ['goagain' num2str(j)];
        disp(folder_name);
        
        %% if not the zeroth folder then read in the final tracked positions of all beads from the previous folder.
        if i == folder_start
            %do nothing
        else
            xinitguess = csvread([directory date '/' ['goagain' num2str(j-1)] '/' 'xlastframe' '.dat']);
            yinitguess = csvread([directory date '/' ['goagain' num2str(j-1)] '/' 'ylastframe' '.dat']);
        end
        
        %% this bit sorts out the not to track ridiculousness
        for g = 0:i % in given folder number
                R = find(wbag == g);
                not_to_track = [not_to_track R];
        end
                
        %%
        [xqi,yqi] = trackerclean_9_0(directory,date,folder_name,descriptor,roi_size,2,no_frames,start,no_steps,frame_int,xinitguess,yinitguess,descriptwo,frame_rate,...
            no_it,...
            no_points_z_porabola,...
            truncate,...
            cut_off,...
            m,...
            frame_ratio,...
            no_beads,...
            startbead,...
            endbead,...
            not_to_track);
        
        xinitguess = xqi(end,:)-(roi_size/2);
        yinitguess = yqi(end,:)-(roi_size/2);
        
        
        csvwrite([directory date '/' folder_name '/' 'xlastframe' '.dat'],xinitguess);
        csvwrite([directory date '/' folder_name '/' 'ylastframe' '.dat'],yinitguess);
        
        
        
    end
end
end

function [xqi,yqi] = trackerclean_9_0(path,date,subfolder,descriptor,roisize,zLUT,no_frames,start,no_steps,frame_int,xer,yer,descriptwo,frame_rate,...
    no_it,...
    no_points_z_porabola,...
    truncate,...
    cut_off,...
    m,...
    frame_ratio,...
    no_beads,...
    startbead,...
    endbead,...
    to_track)
% zLUT = 0 do nothing
% zLUT = 1 measure LUT
% zLUT = 2 read in LUT

% be careful here as to what you are actually tracking - it's bloody confusing because of the
% way the list works. there's no correlation bewtween list order and bead number, the whole p is a member
% of thing is to fix it but it depends if you want to include or exclude.
% If academia valued robust code over publishable results then I would fix it but I don't have time.

disp(frame_ratio);
x = xer;
y = yer;

for p = 1:no_beads
    if ismember(p,to_track) == 1   % if current bead is a member of array 'to_track' then do nothing - i.e. the numbers that are absent are tracked
        x(p) = 250;
        y(p) = 250;       
    else       
        % do nothing
    end
end

x_w_roi = zeros(1,no_beads);
y_w_roi = zeros(1,no_beads);

% deltax_refined = zeros(1,no_beads);
% deltay_refined = zeros(1,no_beads);

xfromqi = zeros(1,no_beads);
yfromqi = zeros(1,no_beads);

xqi = zeros(frame_ratio,no_beads);
yqi = zeros(frame_ratio,no_beads);

LUT = zeros(no_steps,floor(roisize*m - (m-0.5))-truncate-cut_off,no_beads);

hmmmm1 = zeros(frame_ratio,floor(roisize*m - (m-0.5)),no_beads);

radial_profile_new = zeros(frame_ratio,floor(roisize*m - (m-0.5))-truncate-cut_off,no_beads);

LUT_new = zeros(frame_ratio,floor(roisize*m - (m-0.5))-truncate-cut_off,no_beads);

diff = zeros(no_steps,floor(roisize*m - (m-0.5))-truncate-cut_off,no_beads);
fit_quad_arrayz = zeros(no_points_z_porabola,no_beads);
pz = zeros(3,no_beads);

refinedz = zeros(frame_ratio,no_beads);
refinedzinmicrons = zeros(frame_ratio,no_beads);

for frame_no_i = 1:frame_ratio
    %     clear raw_img
    if mod(frame_no_i-1,10) == 0    %print out folder number every now and again
        disp(subfolder);
    else
        %do nothing
    end
    frame_no = ((frame_no_i-1) * frame_int) + start;
    disp(frame_no);
    
    if zLUT == 2
        readinname = [path date '/' subfolder '/' sprintf('frame%1d.jpg', frame_no)];
    elseif zLUT == 1
        readinname = [path date '_LUT/' sprintf('zstack%03d.jpg', frame_no)];
    elseif zLUT ~= 1 && zLUT ~= 2
        readinname = [path date '/' sprintf('frame%1d.jpg', frame_no)];
    end
    
    FOV = double(imread(readinname))/255;
    
    %% if you uncomment this then you can see every frame/FOV dispalyed on screen
    %                 figure(56)
    %                 clf
    %                 image(FOV*255)
    %                 hold on
    %                 colormap(gray(256))
    %                 plot(x,y,'bo')
    %                 drawnow
    
    %%
    for p = startbead:endbead
        %% if you uncomment this then you can see which bead it is currently running through.
%                 bead_is = p
        
        %%
        if ismember(p,to_track) == 1 % if current bead is a member of array 'to_track' then do nothing - i.e. the numbers that are absent are tracked
            xfromqi(p) = (roisize/2) + 1;
            yfromqi(p) = (roisize/2) + 1;
            hmmmm1(frame_no_i,:,p) = 25;  
        else
            
            B_dig = imcrop(FOV, [x(p), y(p), roisize-1, roisize-1]);
            
            %% if you uncomment this then you can see a particular bead (p) displayed on the screen
            %             if p==3
            %                 figure(57)
            %                 image(B_dig*255)
            %                 hold on
            %                 colormap(gray(256))
            %                 drawnow
            %             end
            
            meanpad = mean2(B_dig);
            B_dig_pad = padarray(B_dig,[10 10],meanpad);
            
            [x_w_roi(p),y_w_roi(p)] = cofm(B_dig);
            
            guessx = x_w_roi(p);
            guessy = y_w_roi(p);
            
            for i = 1:no_it
                
                [guessx,guessy,hmmmm1(frame_no_i,:,p)] = qi4(B_dig_pad,guessx,guessy,roisize,m);
                xfromqi(p) = guessx;
                yfromqi(p) = guessy;
                
            end
        end
    end
    
    xqi(frame_no_i,:) = xfromqi+x-1;   %minus one becuse of the weird way imcrop works
    yqi(frame_no_i,:) = yfromqi+y-1;   %minus one becuse of the weird way imcrop works
    
    x = round((xfromqi+x-1)-(roisize/2)); %make roi drift or not
    y = round((yfromqi+y-1)-(roisize/2)); %make roi drift or not
    
end

for g = 1:frame_ratio
    for u = startbead:endbead
        if u==0
        else
            LUT_new(g,:,u) = hmmmm1(g,cut_off+1:end-truncate,u);
            radial_profile_new(g,:,u) = hmmmm1(g,cut_off+1:end-truncate,u);
        end
    end
end

for v = startbead:endbead
    p = v;
    if zLUT == 1
        flippedhmmmm1 = flip(LUT_new,1);
        
        
        filename2 = [path date '_LUT/' '' descriptor 'bead_no' num2str(v) 'new.dat'];
        
        
        dlmwrite(filename2,flippedhmmmm1(:,:,v),'precision','%.6f');
        
        refinedzinmicrons = 0;
        
        refinedz = zeros(no_frames,no_beads);
        filenametosave = [path date '\' '' descriptor 'xyzdatasailor.dat'];
        
    elseif zLUT == 2
        if ismember(p,to_track) == 1
            LUT(:,:,v) = 25;
            radial_profile_new(:,:,v) = 25;
        else
            aaa = [path date '_LUT' '/' descriptor 'bead_no' num2str(v) 'LUT.dat'];
            disp(aaa);
            
            LUT(:,:,v) = csvread([path date '_LUT' '/' descriptor 'bead_no' num2str(v) 'new.dat']);
            
            LUT = bsxfun(@minus, LUT(:,:,v), min(LUT(:,:,v),[],2));
            LUT = bsxfun(@times, LUT, 1./((max(LUT,[],2))-min(LUT,[],2)));
            LUT = bsxfun(@minus, LUT, mean(LUT(:,end-10:end),2));
            
            rad_prof_norm = bsxfun(@minus, radial_profile_new, min(radial_profile_new,[],2));
            rad_prof_norm = bsxfun(@times, rad_prof_norm, 1./((max(rad_prof_norm,[],2))-min(rad_prof_norm,[],2)));
            rad_prof_norm = bsxfun(@minus, rad_prof_norm, mean(rad_prof_norm(:,end-10:end,:),2));
        end
        
        for h = 1:frame_ratio     
            if ismember(p,to_track) == 1   
                refinedz(h,v) = 25;
            else
                
                % LUT is no_frames by profile where frame is step
                % rad_prof_norm is no_frames by profile by no_beads
                % diff is steps by profile by no_beads (in this case for a given frame h)
                % testing is steps by summed difference (single value) by no_beads (in this case for a given frame h)
                
                diff(:,:,v) = bsxfun(@minus, LUT(:,1:end), rad_prof_norm(h,1:end,v));
                
                % testing = sum(abs(diff),2);
                % alternatively you can use
                testing = sum(diff.^2,2);
                [~,gz] = min(testing,[],1);
                
                fit_quad_arrayz(:,v) = linspace(gz(:,:,v)-((no_points_z_porabola-1)/2),gz(:,:,v)+((no_points_z_porabola-1)/2),no_points_z_porabola);
                vectorz = linspace(1,size(LUT,1),size(LUT,1)).';
                
                testy = fit_quad_arrayz(:,v) <= 0 | fit_quad_arrayz(:,v) >= size(LUT,1);
                testys = any(testy);
                
                if testys == 0
                    pz(:,v) = danpolyfit(vectorz(fit_quad_arrayz(:,v)),testing(fit_quad_arrayz(:,v),1,v),2);    
                elseif testys == 1
                    pz(:,v) = NaN;   
                end
                
                refinedz(h,v) = -pz(2,v)./(2.*pz(1,v));
                refinedzinmicrons(h,v) = (refinedz(h,v)-1)*0.2;             
            end
        end  
        filenametosave = [path date '/' subfolder '/' '' descriptor descriptwo 'xyzdata.dat'];
    else
        % do nothing
        filenametosave = [path date '/' subfolder '/' '' descriptor descriptwo 'xyzdata.dat'];
    end
end


a = xqi';
size(a)
b = yqi';
size(b)
c = refinedz.';
size(c)
col_interleave = reshape([a(:) b(:) c(:)]',3*size(a,1), [])';
% bbb = filenametosave
dlmwrite(filenametosave,col_interleave,'newline','pc','precision','%.6f');

figure(95)
subplot(3,1,1)
plot(linspace(1,no_frames,(no_frames)/frame_int)*(1/frame_rate),xqi(:,1),'b-')
xlabel('time (seconds)')
ylabel('x position (pixels)')

subplot(3,1,2)
plot(linspace(1,no_frames,(no_frames)/frame_int)*(1/frame_rate),yqi(:,1),'b-')
xlabel('time (seconds)')
ylabel('y position (pixels)')

subplot(3,1,3)
plot(linspace(1,no_frames,(no_frames)/frame_int)*(1/frame_rate),refinedz(:,1),'b-')
xlabel('time (seconds)')
ylabel('z position (pixels)')
end

function [final_x,final_y,hmmmm1] = qi4(image_in,nextguessx,nextguessy,roisize,m)
%% quadrant interpolation algorithm

firstguessx = nextguessx;
firstguessy = nextguessy;

no_po_r = roisize*m - (m-0.5);
sam = 1/(m*2);
samtestx = linspace(sam, ((m*2)-1)*sam,(m*2)-1);
samtesty = linspace(sam, ((m*2)-1)*sam,(m*2)-1);
interp_coordsx = linspace(1,roisize,no_po_r*2);
interp_coordsy = linspace(1,roisize,no_po_r*2);

no_theta = 20;
klin = linspace(1,no_theta,no_theta);
thetalin = (klin)*(90/(no_theta+1))*pi/180;
thetaminuslin = (klin)*(-90/(no_theta+1))*pi/180;

%%%% change this back to interp_profile (no 2) to return to slower correct
%%%% version
[ql,qr,hmmmm1] = interp_profile5(image_in,firstguessx,firstguessy,0,thetalin,thetaminuslin,interp_coordsx,no_po_r,1,roisize);
[qt,qb,~] = interp_profile5(image_in,firstguessx,firstguessy,1,thetalin,thetaminuslin,interp_coordsy,no_po_r,2,roisize);

if mod(firstguessx,1)~=0 && mod(firstguessx,1)~=0.5 %The point 5 is an ugly fix for a strange bug
    x_avg_profile = [ql(1:end) qr];
elseif mod(firstguessx,1)==0 || mod(firstguessx,1)==0.5
    x_avg_profile = [ql(1:end-1) qr];
end

if mod(firstguessy,1)~=0
    y_avg_profile = [qt(1:end) qb];
else
    y_avg_profile = [qt(1:end-1) qb];
end

testxtest = ones(1,2*(m-0.5))*x_avg_profile(1);
testytest = ones(1,2*(m-0.5))*y_avg_profile(1);

Cx = oneDcrosscorr([testxtest x_avg_profile]);
Cy = oneDcrosscorr([testytest y_avg_profile]);

[~,gx] = max(Cx);
[~,gy] = max(Cy);

if gx <= 2
    gx = 3;
elseif gx >= (2*m*roisize)-1
    gx = (2*m*roisize)-2;
end

if gy <= 2
    gy = 3;
elseif gy >= (2*m*roisize)-1
    gy = (2*m*roisize)-2;
end

fit_quad_arrayx = linspace(gx-2,gx+2,5);
fit_quad_arrayy = linspace(gy-2,gy+2,5);

vectorx = [samtestx interp_coordsx];
vectory = [samtesty interp_coordsy];

[rx,~,mux] = danpolyfit(vectorx(fit_quad_arrayx),Cx(fit_quad_arrayx),2);

Bx = (rx(2)*mux(2) - 2*rx(1)*mux(1))/(mux(2).^2);
Ax = rx(1)/(mux(2).^2);

refinedx = -Bx./(2.*Ax);

delta_rx = (refinedx-(roisize/2))/2;
final_x = (roisize/2) + delta_rx;

[ry,~,muy] = danpolyfit(vectory(fit_quad_arrayy),Cy(fit_quad_arrayy),2);

By = (ry(2)*muy(2) - 2*ry(1)*muy(1))/(muy(2).^2);
Ay = ry(1)/(muy(2).^2);

refinedy = -By./(2.*Ay);

delta_ry = (refinedy-(roisize/2))/2;
final_y = (roisize/2) + delta_ry;

end

function [ql,qr,hmmmm] = interp_profile5(ROI,xguessin,yguessin,trans,thetalin,thetaminuslin,interp_coords,no_po_r,skip,roi_size)
%interpolation of profiles
roi_size = roi_size + 20;

if trans == 0
    xguess = xguessin;
    yguess = yguessin;
elseif trans == 1
    xguess = yguessin;
    yguess = xguessin;
end

if trans == 0
    img = ROI;
elseif trans == 1
    img = ROI.';
end
% sz_a = size(interp_coords)
% sz_b = size(sin((pi/2)-thetalin).')

% x45 = bsxfun(@times,(interp_coords-xguess),sin((pi/2)-thetalin).') + xguess + 10; % the 10 is due to me adding on a boundary that is 10 all the way around
% y45 = bsxfun(@times,(interp_coords-xguess),sin(thetalin).') + yguess + 10;

% sz_c = size(x45)
x45 = (sin((pi/2)-thetalin).' * (interp_coords-xguess)) + xguess + 10;
y45 = (sin(thetalin).' * (interp_coords-xguess)) + yguess + 10;
% sz_d = size(test_x45)

xzero = floor(x45); %+ ceil(xcom)
xone = ceil(x45); %+ ceil(xcom)
yzero = floor(y45); %+ ceil(xcom)
yone = ceil(y45); %+ ceil(xcom)

speedtest = yone-y45;
speedtest2 = y45-yzero;
speedtest3 = xone-x45;
speedtest4 = x45-xzero;

idxzz = yzero + (xzero-1)*roi_size;  %%OMG this is faster than sub2ind
idxzo = yzero + (xone-1)*roi_size;
idxoz = yone + (xzero-1)*roi_size;
idxoo = yone + (xone-1)*roi_size;

f_intereped_total = (img(idxzz).*(speedtest3).*(speedtest) + img(idxzo).*(speedtest4).*(speedtest) + img(idxoz).*(speedtest3).*(speedtest2) + img(idxoo).*(speedtest4).*(speedtest2));

%%%%

% xminus45 =  bsxfun(@times,(interp_coords-xguess),sin((pi/2)-thetaminuslin).') + xguess + 10;
% yminus45 =  bsxfun(@times,(interp_coords-xguess),sin(thetaminuslin).') + yguess + 10;

xminus45 =  (sin((pi/2)-thetaminuslin).' * (interp_coords-xguess)) + xguess + 10;
yminus45 =  (sin(thetaminuslin).' * (interp_coords-xguess)) + yguess + 10;

xzerominus = floor(xminus45); %+ ceil(xcom)
xoneminus = ceil(xminus45); %+ ceil(xcom)
yzerominus = floor(yminus45); %+ ceil(xcom)
yoneminus = ceil(yminus45); %+ ceil(xcom)

speedtestminus = yoneminus-yminus45;
speedtest2minus = yminus45-yzerominus;
speedtest3minus = xoneminus-xminus45;
speedtest4minus = xminus45-xzerominus;

idxzzminus = yzerominus + (xzerominus-1)*roi_size;  %%OMG this is faster than sub2ind
idxzominus = yzerominus + (xoneminus-1)*roi_size;
idxozminus = yoneminus + (xzerominus-1)*roi_size;
idxoominus = yoneminus + (xoneminus-1)*roi_size;

f_intereped_total2 = (img(idxzzminus).*(speedtest3minus).*(speedtestminus) + img(idxzominus).*(speedtest4minus).*(speedtestminus) + img(idxozminus).*(speedtest3minus).*(speedtest2minus) + img(idxoominus).*(speedtest4minus).*(speedtest2minus));

avg_prof = mean(f_intereped_total(2:end-1,:),1);
avg_prof2 = mean(f_intereped_total2(2:end-1,:),1);

argh(1,:) = interp_coords;

I = find(argh(1,:) >= min(xguess));

J = find(argh(1,:) <= min(xguess));

br = avg_prof(I);
tl = avg_prof(J);
tr = avg_prof2(I);
bl = avg_prof2(J);

smax = max([size(br) size(tl) size(tr) size(bl)]);
if smax - size(br,2) == 0
    a = [];
else
    a = ones((smax - size(br,2)),1);
end

if smax - size(tl,2) == 0
    b = [];
else
    b = ones((smax - size(tl,2)),1);
end

if smax - size(tr,2) == 0
    c = [];
else
    
    c = ones((smax - size(tr,2)),1);
end

if smax - size(bl,2) == 0
    d = [];
else
    
    d = ones((smax - size(bl,2)),1);
end

if skip == 2
    hmmmm = 0;
else
    hmmmm = [br a.'*br(end)] + fliplr([b.'*tl(1) tl]) + [tr c.'*tr(end)] + fliplr([d.'*bl(1) bl]);
    hmmmm = hmmmm(1:floor(no_po_r));
end

qr = br + tr;
ql = tl + bl;

end

function [x,y] = cofm(cr_dig)

cr_dig_lin = reshape(cr_dig,[1,numel(cr_dig)]);

D_neg = (cr_dig(:,:)-median(cr_dig_lin));   %subtract median off cropped image from cropped image

thresh = abs(D_neg);      %make all negative values zero

% coljx = zeros(1,size(thresh,1));    %pre-assign for speed
% rowix = zeros(1,size(thresh,2));    %pre-assign for speed
% rowiy = zeros(1,size(thresh,1));    %pre-assign for speed
% coljy = zeros(1,size(thresh,2));    %pre-assign for speed

testj = linspace(1,size(thresh,2),size(thresh,2));
newcoljx = bsxfun(@times,testj,thresh);
newrowix = sum(newcoljx,2);

%calculate com in x

x = sum(newrowix)./sum(sum(thresh));
%
testi = linspace(1,size(thresh,1),size(thresh,1));
newrowiy = bsxfun(@times,testi.',thresh);
newcoljy = sum(newrowiy,1);

%calculate com in y

y = sum(newcoljy)./sum(sum(thresh));

end

function [xcorr1D] = oneDcrosscorr(profile)

xcorr1D = fftshift(ifft(fft(profile).*((fft(fliplr(profile))).'')));

end

function [p,S,mu] = danpolyfit(x,y,n)
%   This is dan's altered polyfit program that removes overheads for speed

%POLYFIT Fit polynomial to data.
%   P = POLYFIT(X,Y,N) finds the coefficients of a polynomial P(X) of
%   degree N that fits the data Y best in a least-squares sense. P is a
%   row vector of length N+1 containing the polynomial coefficients in
%   descending powers, P(1)*X^N + P(2)*X^(N-1) +...+ P(N)*X + P(N+1).
%
%   [P,S] = POLYFIT(X,Y,N) returns the polynomial coefficients P and a
%   structure S for use with POLYVAL to obtain error estimates for
%   predictions.  S contains fields for the triangular factor (R) from a QR
%   decomposition of the Vandermonde matrix of X, the degrees of freedom
%   (df), and the norm of the residuals (normr).  If the data Y are random,
%   an estimate of the covariance matrix of P is (Rinv*Rinv')*normr^2/df,
%   where Rinv is the inverse of R.
%
%   [P,S,MU] = POLYFIT(X,Y,N) finds the coefficients of a polynomial in
%   XHAT = (X-MU(1))/MU(2) where MU(1) = MEAN(X) and MU(2) = STD(X). This
%   centering and scaling transformation improves the numerical properties
%   of both the polynomial and the fitting algorithm.
%
%   Warning messages result if N is >= length(X), if X has repeated, or
%   nearly repeated, points, or if X might need centering and scaling.
%
%   Class support for inputs X,Y:
%      float: double, single
%
%   See also POLY, POLYVAL, ROOTS, LSCOV.

%   Copyright 1984-2011 The MathWorks, Inc.

% The regression problem is formulated in matrix format as:
%
%    y = V*p    or
%
%          3  2
%    y = [x  x  x  1] [p3
%                      p2
%                      p1
%                      p0]
%
% where the vector p contains the coefficients to be found.  For a
% 7th order polynomial, matrix V would be:
%
% V = [x.^7 x.^6 x.^5 x.^4 x.^3 x.^2 x ones(size(x))];

if ~isequal(size(x),size(y))
    error(message('MATLAB:polyfit:XYSizeMismatch'))
end

x = x(:);
y = y(:);

if nargout > 2
    mu = [mean(x); std(x)];
    x = (x - mu(1))/mu(2);
end

% Construct Vandermonde matrix.
V(:,n+1) = ones(length(x),1,class(x));
for j = n:-1:1
    V(:,j) = x.*V(:,j+1);
end

% Solve least squares problem.
[Q,R] = qr(V,0);
% ws = warning('off','all');
p = R\(Q'*y);    % Same as p = V\y;
% % % % warning(ws);
% % % % if size(R,2) > size(R,1)
% % % %    warning(message('MATLAB:polyfit:PolyNotUnique'))
% % % % elseif warnIfLargeConditionNumber(R)
% % % %     if nargout > 2
% % % %         warning(message('MATLAB:polyfit:RepeatedPoints'));
% % % %     else
% % % %         warning(message('MATLAB:polyfit:RepeatedPointsOrRescale'));
% % % %     end
% % % % end

if nargout > 1
    r = y - V*p;
    % S is a structure containing three elements: the triangular factor from a
    % QR decomposition of the Vandermonde matrix, the degrees of freedom and
    % the norm of the residuals.
    S.R = R;
    S.df = max(0,length(y) - (n+1));
    S.normr = norm(r);
end

p = p.';          % Polynomial coefficients are row vectors by convention.

% % % % function flag = warnIfLargeConditionNumber(R)
% % % % if isa(R, 'double')
% % % %     flag = (condest(R) > 1e+10);
% % % % else
% % % %     flag = (condest(R) > 1e+05);
% % % % end
end

function [x,y] = findpositionsofbeads(directory,date,imgname,descriptor,roisize,cc_thresh,rotate)

Z = imread([directory date '/' imgname]);   %read in image
A=imrotate(Z,rotate);

imshow(A);        %display image
h = imrect(gca,[10 10 roisize roisize]);    %draw an rectangle interactive rectangle
setResizable(h,0);                          %prevent the rectangle from being changed in size
pos = wait(h);      %wait for user input
delete(h);          %delete handle
close

xmin = pos(1);      %starting x coord of rectangle
ymin = pos(2);      %strating y coord of rectangle

B = imcrop(A, [xmin, ymin, roisize-1, roisize-1]); %create an image (template) that contains only the rectangle image

A_dig = im2double(A);   %double precision version of master image

%A_dig=A_dig-median(reshape(A_dig,[1,Adig_size(1)*Adig_size(2)]));
B_dig = im2double(B);   %double precision version of template image

% figure(88)
% image(B_dig*255)
% colormap(gray(256))
% axis equal tight

%%%% I have to decide whether to subtract some kind of background
%A_dig=A_dig-median(reshape(A_dig,[1,Adig_size(1)*Adig_size(2)]));
%C = reshape(B_dig,[1,Bdig_size(1)*Bdig_size(2)]);
%B_dig=B_dig-median(C);

%%%% this function only works with the signal processing toolbox
%twoDCC = xcorr2(A_dig,B_dig);

%%%% this fucntion works with the image processing toolbox
twoDCC = normxcorr2(B_dig,A_dig);   %2D cross correlation
scc=size(twoDCC);                   %how big is it?
[max_cc, ~] = max(max(twoDCC));   %what's the maximum value in it?

twoDCC=(twoDCC)/max_cc;         %normalise

indices = find(twoDCC > cc_thresh);     %find pixel positions where the cross correlation is larger than a certain threshold

twoDCC(indices) = 1;            %set the pixel values given by the elements in indices to unity
twoDCC(find(twoDCC<1))=0;       %set everything else to zero

roi_half = round(roisize/2);

%this for loop is a dirty way of turning blobs into single pixel points
for i = 2:size(twoDCC,1)
    for j = 2:size(twoDCC,2)
        if twoDCC(i,j) == 1
            twoDCC((-roi_half:roi_half)+i,(-roi_half:roi_half)+j) = 0;
            twoDCC(i,j) = 1;
        end
    end
end

indices_one = find(twoDCC == 1);    %find all the indices of where the points in the 2D xcor are
[U,V] = ind2sub(scc,indices_one);   %get array subscripts from indices

shift_x = (size(B_dig,1) - 1)/2;     %xcor creates a larger image so it has to be compensated for
shift_y = (size(B_dig,2) - 1)/2;

U = U - shift_y;
V = V - shift_x;

x_w_roi = zeros(1,length(U));        %pre-assign for speed
y_w_roi = zeros(1,length(U));        %pre-assign for speed

for b = 1:length(U)
    template_n = imcrop(A, [V(b)-roi_half, U(b)-roi_half, roisize - 1, roisize - 1]);
    template_n_dig = im2double(template_n);
%     figure(5555)
%     image(template_n_dig*255)
%     colormap(gray(256))
    [x_w_roi(b),y_w_roi(b)] = cofm(template_n_dig);
end

figure(8)
image(A_dig*255)
colormap(gray(256))
axis equal tight

xfromcom = round(V-roi_half+(x_w_roi.'-roi_half));
yfromcom = round(U-roi_half+(y_w_roi.'-roi_half));

for b = 1:length(U)
    template_n2 = imcrop(A, [xfromcom(b), yfromcom(b), roisize - 1, roisize - 1]);
    template_n2_dig = im2double(template_n2);
    deltax_refined(b) = subpixelfit(template_n2_dig,roisize);
    deltay_refined(b) = subpixelfit(template_n2_dig.',roisize);
end
size(xfromcom);
size(deltax_refined);
x = round(deltax_refined + xfromcom.' - 1);
y = round(deltay_refined + yfromcom.' - 1);

for q = 1:length(U)
    rectangle('Position', [x(q)-roi_half, y(q)-roi_half, roisize, roisize], 'EdgeColor', [1 1 1]);
end

csvwrite([directory date '/' 'csvlistx' descriptor '.dat'],x)
csvwrite([directory date '/' 'csvlisty' descriptor '.dat'],y)

figure(56)
image(A_dig*255)
hold on
colormap(gray(256))
plot(V,U,'wo')
hold on
plot(x,y,'bo')

end

function [refined_pos] = subpixelfit(B_dig,roisize)

P_avg = mean_prof(B_dig,roisize);
C = oneDcrosscorr(P_avg);


[~,g] = max(C);

if g <= 2
    g = 3;
elseif g >= roisize-1
    g = roisize-2;
end

fit_quad_array = linspace(g-2,g+2,5);
p = danpolyfit(fit_quad_array,C(fit_quad_array),2);
refined = -p(2)./(2.*p(1));
refined_pos = (roisize/2)+((refined-(roisize/2))/2);

end

function [P_avg] = mean_prof(cropped_image,roisize)

bottom = round(0.2*roisize);
top = round(0.8*roisize);
P = zeros(top-bottom+1,size(cropped_image,2));

t = linspace(bottom,top,top-bottom+1);
ind = t-bottom+1;
P(1:length(ind),:) = cropped_image(t,:);

P_avg = mean(P,1);
end

function displaybeads(directory,date,imgname,descriptor,roisize)

x_auto = csvread([directory date '/' 'csvlistx' descriptor '.dat']);
y_auto = csvread([directory date '/' 'csvlisty' descriptor '.dat']);
Z = imread([directory date '/' imgname]);    %read in image

A=Z;

A_dig = im2double(A);   %double precision version of master image

roi_half = round(roisize/2);

figure(8)
image(A_dig*255)
colormap(gray(256))
axis equal tight

for q = 1:length(x_auto)
    rectangle('Position', [x_auto(q)-roi_half, y_auto(q)-roi_half, roisize, roisize], 'EdgeColor', [1 1 1]);
    text(x_auto(q)+roi_half, y_auto(q)+roi_half, num2str(q), 'color', 'white')
end

end

function [x_add,y_add] = addpos(directory,date,imgname,descriptor,roisize)
[directory date '/' 'csvlistx' descriptor '.dat'];
x_auto = csvread([directory date '/' 'csvlistx' descriptor '.dat']);
y_auto = csvread([directory date '/' 'csvlisty' descriptor '.dat']);

% single jpeg
Z = imread([directory date '/' imgname]);    %read in image

%A = rgb2gray(Z);        %if needed make the image greyscale
A=Z;

imshow(A);        %display image

for q = 1:length(x_auto)
    rectangle('Position', [x_auto(q)-(roisize/2), y_auto(q)-(roisize/2), roisize, roisize], 'EdgeColor', [1 1 1]);
end

h = imrect(gca,[10 10 roisize roisize]);    %draw an rectangle interactive rectangle
setResizable(h,0);                          %prevent the rectangle from being changed in size
pos = wait(h);      %wait for user input
delete(h);          %delete handle

close
roi_half = round(roisize/2);
xmin = pos(1)      %starting x coord of rectangle
ymin = pos(2)      %strating y coord of rectangle

A_dig = im2double(A);   %double precision version of master image



template_n = imcrop(A, [xmin, ymin, roisize-1, roisize-1]);
template_n_dig = im2double(template_n);
[x_w_roi,y_w_roi] = cofm(template_n_dig)

figure(8)
image(A_dig*255)
colormap(gray(256))
axis equal tight

x_add=round(xmin+(x_w_roi.'));
y_add=round(ymin+(y_w_roi.'));

x_final = [x_auto floor(x_add)];
y_final = [y_auto floor(y_add)];

for q = 1:length(x_final)
    rectangle('Position', [x_final(q)-(roisize/2), y_final(q)-(roisize/2), roisize, roisize], 'EdgeColor', [1 1 1]);
end

csvwrite([directory date '\' 'csvlistx' descriptor '.dat'],x_final);
csvwrite([directory date '\' 'csvlisty' descriptor '.dat'],y_final);

figure(56)
image(A_dig*255)
hold on
colormap(gray(256))
hold on
plot(x_add,y_add,'bo')

end