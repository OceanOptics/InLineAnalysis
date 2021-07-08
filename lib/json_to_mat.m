function matfile = json_to_mat(jsonfile)
% Convert InLineAnalaysis user input .json files into new .mat file format 
% Guillaume Bourdin
%
%%
jsonfile = loadjson(jsonfile);
fieldna = fieldnames(jsonfile);
for i = 1:size(fieldna,1)
  if size(jsonfile.(fieldna{i}),2) == 1
    matfile.(fieldna{i}) = datenum(jsonfile.(fieldna{i}){1});
  elseif size(jsonfile.(fieldna{i}),2) == 2
    matfile.(fieldna{i}) = [datenum(jsonfile.(fieldna{i}){1}) datenum(jsonfile.(fieldna{i}){2})];
  else
    error('json file size not supported')
  end
end