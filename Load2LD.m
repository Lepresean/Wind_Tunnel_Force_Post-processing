function [data] = Load2LD(data,theta,chord,span,biasCorr,tunnelCorr)
%this function takes a single data set from raw labview output up through calculating aero coefficients corrected for tunnel effects and load cell bias
%   basically ported over from the force_post_process.m script as it existed 3_10_20
%   theta is angle in deg CCW mini45 x-axis to airfoil chord?
%% remove unnecessary zero rows, check for drift, flag sweep direction, reorder
[data,up_down] = removeZeros(data);
%% apply correction for load cell bias?
if biasCorr == 1
    data = biasCorrection(data);
end
%should Fx,Fy error est be maintained and converted to use as bias estimates?
%% coordinate transform load cell axis to wind axis
data = coordTransform(data,theta); %could change to allow for speed factor, aoa_offset
% data_wind_nbc = data_wind; %save data with no horizontal buoyancy correction
%% horizontal buoyancy correction
if tunnelCorr == 1
    data = horoBuoy(data,chord); %there are a few constants hard coded in here that should be changed for different airfoils/speeds/tunnels
end
%% non-dimensionalize data/calculate aero coeffs
data = nonDim(data,chord,span);
% data_coef_nbc = nonDim(data_wind_nbc,chord,span);
%% calculate and apply wall corrections
if tunnelCorr == 1
    data = wallCorr(data,chord);
end
%% add in sweep direction marker
[row,col]=size(data);
data(1,col+1) = up_down;
end
function [data,up_down] = removeZeros(data)
%find zero rows, checks the drift was not too large, removes zero rows
[row, ~] = size(data);
%% check for drift, remove rows of zero speed
zero_index = 0;
for i=1:row
    if data(i,1) == 0 %if speed is 0
        if abs(data(i,4))> .2 || abs(data(i,5))> .2 || abs(data(i,6))> 1     %if Fx Fy or Fz are greater than assigned limits
            error(strcat('This force data may have drifted too much at 0 speed. Check line: ',num2str(i)));
        else %want to delete these zero rows
            if zero_index == 0
                zero_index = i;
            else
                zero_index=cat(1,zero_index,i); %update marker to indicate that the current row has only zero data
            end
        end
    else
    end
end
if ~zero_index == 0
    data(zero_index,:)=[];
    [row, ~] = size(data);
end
%% check direction of aoa sweep, reorder vectors
%up_down=1 means run sweeps upward
%up_down=0 means run sweeps downward
%up_down=2 means run doesn't sweep upward or downward
if data(2,3)<data(3,3) %if aoa 2 is less than aoa 3
    if data(5,3)<data(6,3) %if aoa 5 < aoa 6 (this should be more robust, but probably fine for now I guess)
        up_down = 1; %flag =1 indicates sweep from low aoa to high
        data(:,:) = flip(data(:,:)); %flips data from this run to match upward sweep data
    else
        up_down = 2; %should not happen unless a mistake was made in data aquisition
    end
elseif data(2,3)>data(3,3) %if aoa 2 is > aoa 3
    if data(5,3)>data(6,3) %if aoa 5 < aoa 6
        up_down = 0; %flag =0 indicates sweep from high aoa to low
    else
        up_down = 2; %should not happen unless a mistake was made in data aquisition
    end
else
    up_down= 2; %should not happen unless a mistake was made in data aquisition
end
%% Check for drift, remove rows of zero angle of attack
%need to add section to check drift for zero aoa rows, average? or throw out outside zeros?
zero_aoa_index = 0;
for i=1:row
    if data(i,3) == 0 %if aoa is 0
        if zero_aoa_index == 0
            zero_aoa_index = i;
            data_0_aoa = data(i,4:9);
        else
            zero_aoa_index=cat(1,zero_aoa_index,i); %update marker to indicate that the current row has only zero data
            data_0_aoa = [data_0_aoa;data(i,4:9)];
        end
    end
end

data_0_aoa_mean = mean(data_0_aoa);
limits = [.2;.2;1;.01;.01;.01]; %drift limits for Fx,Fy,Fz,Tx,Ty,Tz (N, Nmm)

[x,y] = size(data_0_aoa);
for n=1:x %check for drift
    for j=1:y
        if abs(data_0_aoa(n,j)-data_0_aoa_mean(j)) > limits(j) %if Fx Fy or Fz drift more than assigned limits
            error(strcat('This force data may have drifted too much at AOA. Check line: ',num2str(zero_aoa_index(n)),' col: ',num2str(j+3)));
        end
    end
end

if x >= 2 %% make sure there is at least one row to delete
    if x == 3
        if zero_aoa_index(1) == 1 && zero_aoa_index(x) == row
            data(row,:)=[];
            data(1,:)=[];
        else
            data(zero_aoa_index(x),:)=[];
            data(z,:)=[];
            warning('Zero AOA rows exist in interior data');
        end
    elseif x > 3
        data(zero_aoa_index(x),:)=[];
        data(z,:)=[];
        num_to_delete = x-3; %number of rows to get rid of (zero aoa rows minus 2 for front and back minus 1 to keep)
        for i= num_to_delete : -1 : 0
            data(zero_aoa_index(x-i),:)=[];
        end
        warning('Zero AOA rows exist in interior data, more than 3 zero aoa rows');
    elseif x == 2
        if zero_aoa_index(x) == row
            data(zero_aoa_index(x),:)=[];
        elseif zero_aoa_index(1) == 1
            data(1,:)=[];
        else
            data(zero_aoa_index(x),:)=[];
        end
    end
end
end
function [data_fix] = biasCorrection(data)
%Calculates the load cell bias in a data set
[row,col]=size(data);
res_angle = atan2d(data(:,5),data(:,4)); %angle CCW from positive x axis on load cell to resultant force measured
% res_angle = reshape(res_angle,row,1);
res = sqrt(data(:,4).^2+data(:,5).^2); %resultant force
% res = reshape(res,row,1);
data_fix=data;
data_fix(:,4) = zeros(row,1);
data_fix(:,5) = zeros(row,1);

[Diff_Fx_est,Diff_Fy_est] = Mini45_Fx_Fy_error(res,res_angle);
data_fix(:,4) = data(:,4)+Diff_Fx_est;
data_fix(:,5) = data(:,5)+Diff_Fy_est;
% data = cat(1,data,data_fix); %add corrected data to end of raw data to do all processing at once
end
function [data_wind] = coordTransform(data,theta)
%Performs coordinate transform on load cell data to convert to wind axis (Lift and Drag)
%   theta is angle in deg CCW mini45 x-axis to airfoil chord?

%AOA + theta = transform_angle between load cell x and wind x
%mini45 -y is aligned to airfoil x (airfoil x out LE)
%mini40 x is 30 CW to airfoil x
%theta defined in initialization section

[row,~]=size(data);

%column order: speed, AOA, Drag, Lift, Mom_z
data_wind = zeros(row,6);   %initialize new matrix for wind axis values
data_wind(:,1) = data(:,1);%*speed_factor; %copy speed column, apply correction to data from 11_17_19
data_wind(:,2) = data(:,2);% copy density column
data_wind(:,3) = data(:,3);%+aoa_offset; %copy aoa column
for i=1:row
    %         aoa = data(i,3)+aoa_offset; %angle between airfoil x and wind x
    %         data_wind(i,3,k) = data(i,4,k)*sind(aoa)-data(i,5,k)*cosd(aoa); %drag = Fx*sin(aoa)-Fy*cos(aoa) Mini45 SAME as formula 3 when theta=90
    %         data_wind(i,4,k) = data(i,4,k)*cosd(aoa)+data(i,5,k)*sind(aoa); %lift = Fx*cos(aoa)+Fy*sin(aoa) Mini45 SAME as formula 4 when theta=90
    data_wind(i,4) = data(i,4)*cosd(theta-data_wind(i,3))-data(i,5)*sind(theta-data_wind(i,3)); %drag, formula 3
    data_wind(i,5) = data(i,4)*sind(theta-data_wind(i,3))+data(i,5)*cosd(theta-data_wind(i,3)); %lift, formula 4
    data_wind(i,6) = data(i,9); %moment z
end
end
function [data] = horoBuoy(data,chord)
%applies a correction to L, D due to horozontal buoyancy in the wind tunnel
% see p. 353 "Low Speed Wind Tunnel Testing" by Rae, Pope
dp_dl = -13.8;              %streamwise static pressure gradient (Pa/m) at model location (got this from pressure gradient measured 9.20.19, originally estimated -5 from old data)
Lambda = 0.24;              %airfoil geometry factor for NACA 0012 (graph p. 352)
height = 36*2.54/100;       %wind tunnel "height" in m (height because the airfoil is on its side, really width for us)
sigma = pi^2/48*(chord/height)^2;
Drag_delta = -6*height^2/pi*Lambda*sigma*dp_dl;     %estimated drag effect due to horizontal buoyancy
data(:,4,:)=data(:,4,:)-Drag_delta;       %correct drag measurements for horizontal buoyancy
end
function [data_coef] = nonDim(data,chord,span)
%non-dimensionalize data/calculate aerodynamic coeffs
[row,~]=size(data);
rho = mean(data(:,2));                 %average density for each run
q = .5*rho*mean(data(:,1))^2;  %average dynamic pressure for each run

data_coef = zeros(row,5);   %initializing array for non-dim data
data_coef(:,1:3) = data(:,1:3);    %copying speed, density, aoa columns
for i=1:row
    data_coef(i,4) = nondim_force(data(i,4),q,chord,span); %section drag coef
    data_coef(i,5) = nondim_force(data(i,5),q,chord,span); %section lift coef
    data_coef(i,6) = nondim_mom(data(i,6),q,chord,span); %section pitching mom coef
end
% column order is now: speed, density, aoa, c_d, c_l, c_m1/4
end
function [data_coef_corr] = wallCorr(data_coef,chord)
% calculate corrections to aerodynamic coefficients for solid blockage, wake blockage, streamline curvature, wake gradient
[row,col]=size(data_coef);
%% solid blockage correction (p. 355 "Low Speed Wind Tunnel Testing" by Rae and Pope)
height = 36*2.54/100;       %wind tunnel "height" in m (height because the airfoil is on its side, really width for us)
sigma = pi^2/48*(chord/height)^2;
volume_foil = 0.7*.12*chord^2*12*2.54/100; %approximate volume of airfoil model (m^3)
K = 0.52;   %value for an airfoil spanning tunnel height (p. 355)
del = .55/100;          %Boundary layer thickness (m) at test location (panel 1) as measured by J. Morris in Fall 2016
del_star = del/6;       %Approx displacement thickness (m) at test location
% A_section1 = (36*2.54/100)*(12*2.54/100); %cross sectional area of tunnel (m^2)
A_section = (36*2.54/100-2*del_star)*(12*2.54/100-2*del_star); %cross sectional area of tunnel minus BL (m^2)
epsilon_sb = K*volume_foil/A_section^(3/2);
%% wake blockage correction and sreamline curvature correction (p. 355, p.359 "Low Speed Wind Tunnel Testing" by Rae and Pope)
epsilon_wb = zeros(row,1);
epsilon = zeros(row,1);
data_coef_corr = zeros(row,col); %new array to store data after solid blockage, wake blockage, and streamline curvature corrections. Consider overwrting data_coef, would that mess up the wake gradient cdo correction?
for i=1:row
    epsilon_wb(i) = chord/height/2*data_coef(i,4);              %wake blockage correction
    epsilon(i) = epsilon_sb+epsilon_wb(i);
    data_coef_corr(i,1) = data_coef(i,1)*(1+epsilon(i));      %correct velocity, note that q changes as well
    data_coef_corr(i,2) = data_coef(i,2); %no correction for density as far as I am aware
    data_coef_corr(i,3) = data_coef(i,3) + 57.3*sigma/2/pi*(data_coef(i,5) + 4*data_coef(i,6));    %streamline curvature correction to aoa (Eqn. 6.20)
    %%% need to check the above formula, was it originally wrong in data_post_process script?
    data_coef_corr(i,5) = data_coef(i,5)*(1-sigma-2*epsilon(i));      %correction to lift coefficient (Eqn. 6.21)
    data_coef_corr(i,6) = data_coef(i,6)*(1-2*epsilon(i))+sigma*data_coef_corr(i,5)/4;  %correction to pitching moment coefficient (Eqn. 6.22)
end
%% separate correction for c_do from the dynamic pressure effect (wake gradient)
%along with the wake gradient term. See Wind Tunnel Testing by Rae, Pope, p. 360, 357
loc_aoa_0 = find(~data_coef(:,3)); %finds the rows that have 0 aoa
[y,~] = size(loc_aoa_0);
if y >= 2 %just to check that nothing weird is being ignored
    error(strcat('Yo, there are more than two zero angle of attack rows in file: ',num2str(k)))
end

c_do = data_coef(loc_aoa_0(1),4);

% if up_down == 1 %if upward sweep
%     %take first zero aoa as c_do position
%     c_do = data_coef(loc_aoa_0(1),3);     %zero lift drag coef is at the first zero aoa row, 3rd col, kth run in data_coef matrix
% elseif up_down == 0 %if downward sweep
%     %take second zero aoa as c_do position
%     c_do = data_coef(loc_aoa_0(2),3);
% else
%     c_do = data_coef(loc_aoa_0(1),3); %just use first time at 0 aoa to define cdo if sweep is wonked
% end

for i=1:row     %apply the c_do correction to all drag coef values
    delta_cdo = -c_do*(3*epsilon_sb+2*epsilon_wb(i)); %zero-lift drag coefficient adjustment for dyn pressure and wake gradient effects
    data_coef_corr(i,4) = data_coef(i,4) + delta_cdo; %drag coefficient adjustment
end
end
function [coef] = nondim_force(force,q,chord,span)
%Turns force into section coefficient
%   D -> c_d
%force = raw force
%chord = airfoil chord
%span = airfoil span
%q = dynamic pressure
coef = force/(q*chord*span);
end
function [coef] = nondim_mom(moment,q,chord,span)
%Turns moment into section coefficient
%   M_z -> c_m
%moment = some moment
%chord = airfoil chord
%q = dynamic pressure
coef = moment/(q*span*chord^2);
end
function [Delta_Fx_est,Delta_Fy_est] = Mini45_Fx_Fy_error(resultant_load,resultant_angle)
%interpolates a table of measured error values for the ATI mini45 load cell to estimate error at given resultant force angle and component load magnitude
%constants copied from load_cell_calibration script, generated from Trial_2 of recalibration campaign 1/30/20
%perhaps a better way would be to create files with these values that this function can reference so it would be capable of dynamic updates, but at this moment I am not expecting these values to need to change much

Fx_error = [-0.0812776024109794,-0.214165846689139,0.210860342211221,0.424754408273175,0.477061255410979,-0.00140301931085995,-0.187919822211220,0.0777273317268250;...
    0.377457340767062,-0.0142375200674167,0.911212619633663,1.35806078881953,1.45908005323294,-0.0314495569325799,-0.654033749633659,-0.241293125819526;...
    1.17026098012314,0.390553708176028,2.00737910247855,3.52760880891223,3.42492622887686,-0.204017682176020,-1.46097191147854,-0.140448239912224];

Fy_error = [0.213803341788779,0.252667669726824,0.0695737395890204,-0.199689249689138,-0.558034265788778,0.0136386002731754,0.380097085410981,0.206904173689139;...
    0.275582448366337,0.568067135180474,0.707057275767062,0.0208954139325863,-1.66158443636634,-0.144154242180472,1.26779501223294,1.09288332106742;...
    0.424174603521447,1.34120406708777,2.92941388712314,1.16660785617604,-3.31349558852145,-0.342914483087775,3.44494000387686,2.35557001882397];

resultant_loads_grid = [4.37068133480380,4.37068133480380,4.37068133480380,4.37068133480380,4.37068133480380,4.37068133480380,4.37068133480380,4.37068133480380;...
    13.1069565225514,13.1069565225514,13.1069565225514,13.1069565225514,13.1069565225514,13.1069565225514,13.1069565225514,13.1069565225514;...
    30.7028085222075,30.7028085222075,30.7028085222075,30.7028085222075,30.7028085222075,30.7028085222075,30.7028085222075,30.7028085222075];

resultant_angles_grid = [-144.017574552423,-98.6731897377625,-59.2927127512379,-10.7725947149021,44.4923767226930,80.1161208145459,126.406311910758,174.110084237203;...
    -144.017574552423,-98.6731897377625,-59.2927127512379,-10.7725947149021,44.4923767226930,80.1161208145459,126.406311910758,174.110084237203;...
    -144.017574552423,-98.6731897377625,-59.2927127512379,-10.7725947149021,44.4923767226930,80.1161208145459,126.406311910758,174.110084237203];

Delta_Fx_gridded = griddedInterpolant(resultant_loads_grid,resultant_angles_grid,Fx_error);
Delta_Fx_est = Delta_Fx_gridded(resultant_load,resultant_angle);

Delta_Fy_gridded = griddedInterpolant(resultant_loads_grid,resultant_angles_grid,Fy_error);
Delta_Fy_est = Delta_Fy_gridded(resultant_load,resultant_angle);
end