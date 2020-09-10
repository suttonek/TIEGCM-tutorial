function [data,units] = get_netcdf_variables(filename,reqvars)


ncid = netcdf.open(filename, 'NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for varid = 0:nvars-1,

	[varname,xtype,dimids,natts] = netcdf.inqVar(ncid,varid);
	if (nargin == 1) || any(strcmp(varname,reqvars))
		eval(sprintf('data.%s = netcdf.getVar(ncid,varid);',varname));
		try eval(sprintf('units.%s = netcdf.getAtt(ncid,varid,''units'');',varname));end
	end

end

%clear ndims nvars ngatts unlimdimid varname xtype dimids natts
netcdf.close(ncid);
