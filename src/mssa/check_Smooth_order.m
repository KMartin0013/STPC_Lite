function [method_order] = check_Smooth_order (Method_Smooth);

%%

Default_Method = ["None","Gaussian300km","Gaussian500km","DDK3","DDK4",...
    "DDK5","DDK6","DDK7"];

Use_method_num = numel(Method_Smooth);

method_order=zeros(1,Use_method_num);
for i=1:Use_method_num

    tf = strcmp(Method_Smooth(i),Default_Method);

    if sum(tf)==0

        disp(['The required filter ' char(Method_Smooth(i)) ' is not supported.'])
        error('Invalid smooth methods.')

    else

        method_order(i)=find(tf==1);

    end
end

end

