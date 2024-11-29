function Tq=CapacityDynamic(k,Delta,T,nTrainHour)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of the capacity on the arc as a function of the time of day
% Equation (30) of the paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FractionCapacity=0.15; % 15 % (percetage) de capacity for freight at peak hours

T24=mod(T,24); % Time of Day

if 10 <= T24 & T24 <=18

    Tq=Delta / (k*nTrainHour *FractionCapacity);

elseif (7 <= T24 & T24 <=10) |  (18 <= T24 & T24 <=24)

    Tq=Delta / (2*k*nTrainHour *FractionCapacity);

else

    Tq=Delta / (k*nTrainHour);

end