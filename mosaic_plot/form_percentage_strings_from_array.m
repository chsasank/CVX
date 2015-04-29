function strs=form_percentage_strings_from_array(array,dec_place)

if ~exist('dec_place','var')
    dec_place=1;
end

format=sprintf('%%1.%df\\n',dec_place);

array=array(:);
array=array/sum(array);
strs=num2str(array*100,format);
strs=[strs ones(length(array),1)*'%'];

strs=cellstr(strs);


