function [labels,rigidBodySegments] = GetJointAndSegmentLabels(templateFilename)

segmentNames = {};
% templateFile = [directory 'modelTemplate.mdh'];
fid = fopen(templateFilename);
block_size = 1;
formatString = '%s %s %s';

index = 1;
while ~feof(fid)
    segarray = textscan(fid, formatString, block_size);
    
    strings = segarray;
    
    if (length(strings) == 3)
        if (isempty(strings{1}) ||  ...
                isempty(strings{2}) || ...
                isempty(strings{3}))
            continue;
        end
        
        if (strcmp(strings{1}{:}, '!') && strcmp(strings{2}{:}, 'Segment') && ~strcmp(strings{3}{:}, 'Info'))
            segmentNames{index} = strings{3}{:};
            index = index + 1;
        end
    end
end
fclose(fid);

fid2 = fopen(templateFilename);
% frewind(fid)

formatString = '%s';
index = 1;
while ~feof(fid2)
    segarray = textscan(fid2, formatString, block_size);
    strings = segarray;
    
    if (length(strings) == 1)
        if (isempty(strings{1}))
            continue;
        end
        
        [ind] = findstr(strings{1}{:}, '/TRACKING_NAMES=');
        
%         fprintf(1,[strings{1}{1} '\n']);
        
        if length(strings{1}{1})>=23
            if (strcmp(strings{1}{1}(1:23),'/REFERENCE_OBJECT_NAMES'))
                s = strrep(strings{1}{1}, '+++', '+');
                s = strrep(s, '++', '+');
                eq = regexp(s,'=');
                pluses = regexp(s,'+');
                starts = [eq+1 pluses+1];
                ends = [pluses-1 length(s)];
                
                for i = 1:length(starts)
                    if length(starts(i):ends(i)) <= 4
                        joints.(segmentNames{index}){i} = s(starts(i):ends(i));
                    end
                end
            end
        end
        
        if (~isempty(ind))
            %       [ind] = findstr(strings{1}{:}, segmentString);
            markerNames = strings{1}{:};
            
            markerNames(1:16) = [];
            markerNameList = regexp(markerNames, '+', 'split');
            
            
            count = 0;
            for k = 1:length(markerNameList)
               if length(markerNameList{k})<=4
                   count = count+1;
                    rigidBodySegments.(segmentNames{index}){count} = markerNameList{k};
               end
            end
            index = index + 1;
        end
        
    end
end

fclose(fid2);


%% Each segment and each joint that shares markers to the adjacent segments has a label
mnames = fieldnames(rigidBodySegments);
for i = 1:length(mnames)
    if strcmp(mnames{i},'RFT')
        label = 1;
    elseif strcmp(mnames{i},'RSK')
        label = 2;
    elseif strcmp(mnames{i},'RTH')
        label = 3;
    elseif strcmp(mnames{i},'LFT')
        label = 4;
    elseif strcmp(mnames{i},'LSK')
        label = 5;
    elseif strcmp(mnames{i},'LTH')
        label = 6;
    elseif strcmp(mnames{i},'RPV')
        label = 7;
    else
        error('Unknown Segment Type %s',mnames{i});
    end
    
    for j = 1:length(rigidBodySegments.(mnames{i}))
        labels.(rigidBodySegments.(mnames{i}){j}) = label;
    end
    
end

%Right ankle markers
for i = 1:length(joints.RFT)
    marker = joints.RFT{i};
   if any(strcmp(marker,joints.RSK))
       labels.(marker) = 8;
   end
end

%Right knee markers
for i = 1:length(joints.RSK)
    marker = joints.RSK{i};
   if any(strcmp(marker,joints.RTH))
       labels.(marker) = 9;
   end
end

%Left ankle markers
for i = 1:length(joints.LFT)
    marker = joints.LFT{i};
   if any(strcmp(marker,joints.LSK))
       labels.(marker) = 10;
   end
end

%Left knee markers
for i = 1:length(joints.LSK)
    marker = joints.LSK{i};
   if any(strcmp(marker,joints.LTH))
       labels.(marker) = 11;
   end
end

