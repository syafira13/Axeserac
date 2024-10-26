function omega = disk_om(r,proret,S)
dr = 1e-4;

dgpp = @(r) (S.g_uv_comp([0,r+dr,pi/2,0],4,4)-S.g_uv_comp([0,r-dr,pi/2,0],4,4))/(2*dr);
dgtt = @(r) (S.g_uv_comp([0,r+dr,pi/2,0],1,1)-S.g_uv_comp([0,r-dr,pi/2,0],1,1))/(2*dr);
dgtp = @(r) (S.g_uv_comp([0,r+dr,pi/2,0],1,4)-S.g_uv_comp([0,r-dr,pi/2,0],1,4))/(2*dr);

omega = (-dgtp(r) + proret*sqrt(dgtp(r)^2 - dgpp(r)*dgtt(r)))/(dgpp(r));