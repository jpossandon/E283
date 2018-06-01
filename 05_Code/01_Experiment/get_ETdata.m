function [data,type] = get_ETdata

datatypes = [3:9,24,25,28,200];
while 1
    type = Eyelink('GetNextDataType');
    if find(datatypes==type)
        data = Eyelink( 'GetFloatData',type);
        break
    end
end
