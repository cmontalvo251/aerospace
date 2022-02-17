function dataOUT = QUANTsim(dataIN,flag,BITS,Dynamic_Range)
%%%%this takes a vector with a certain floating point or fixed
%precision number and degrades it to a precision similar to something found on a
%small mircroprocessor

%%flag = 0 for floating point and 1 for fixed, 2 for dynamic range
%%Floating Point
   %BITS is a vector with structure BITS = [ABITS EBITS]
   %ABITS is the number of bits retained for the mantissa
   %EBITS is the number of bits retained for the exponent integer
%%Fixed Point
   %BITS is a vector with structure BITS = [ABITS EBITS]
   %ABITS is the number of bits retained before the decimal
   %EBITS is the number of bits retained after the decimal
%%Dynamic Range
   %BITS is a scalar specifying the number of partitions to create
   %Dynamic_Range is a vector specifying the max and min of the
   %range of the sensor. If the data point is out of range the max
   %or min value is used. Thus the range of the sensor is split
   %into 2^BITS partitions
dataOUT = dataIN;
if flag == 2
   %%Preliminary Dynamic Range calculations
   upperlimit = Dynamic_Range(2);
   lowerlimit = Dynamic_Range(1);
   quantizedrange = linspace(lowerlimit,upperlimit,2^BITS);
   partitionsize = quantizedrange(2)-quantizedrange(1);
end
[r,c] = size(dataIN);
for i = 1:c
  if flag == 0
    %%Floating Point
    ABITS = BITS(1);
    EBITS = BITS(2);
    for jj = 1:r
       asox = sign(dataIN(jj,i));
       aquant = 2.^mod(log(abs(dataIN(jj,i)))./log(2),1);
       aquant = floor(aquant.*(2.^ABITS))./(2.^ABITS);
       equant = floor(log(abs(dataIN(jj,i)))./log(2));
       esox = sign(equant);
       nquant = floor(log(abs(equant))./log(2)) + 1;
       equant = (2.^(nquant-EBITS)).*floor(equant./(2.^(nquant-EBITS)));
       dataOUT(jj,i) = asox.*aquant.*(2.^(esox.*equant));
    end
  elseif flag == 1
    %%Fixed Point
    for jj = 1:r
       ABITS = BITS(1);
       EBITS = BITS(2);
       asox = sign(dataIN(jj,i));
       aquant = abs(dataIN(jj,i));
       nquant = floor(log(aquant)./log(2)) + 1;
       if ABITS >= nquant
          %%in range
          aquant = floor(aquant.*(2.^EBITS))./(2.^EBITS);
          dataOUT(jj,i) = asox.*aquant;
       else
          %%out of range
          dataOUT(jj,i) = asox.*((2.^ABITS-2.^(-EBITS)));
       end
    end
  elseif flag == 2
     %%Dynamic Range
     for jj = 1:r
        %%Check if out of range
        var = dataIN(jj,i);
        if dataIN(jj,i) >= upperlimit
           dataOUT(jj,i) = upperlimit;
        elseif dataIN(jj,i) <= lowerlimit
           dataOUT(jj,i) = lowerlimit;
        else
           %%In Range
           loc = find(dataIN(jj,i)<quantizedrange,1);
           if dataIN(jj,i) >= quantizedrange(loc-1)+partitionsize/2;
              dataOUT(jj,i) = quantizedrange(loc);
           else
              dataOUT(jj,i) = quantizedrange(loc-1);
           end
        end
     end
  end
end
%disp('Data Quantized')




