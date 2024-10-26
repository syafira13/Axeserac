function risco = disk_r_isco(r_guess,proret,sptm)

dr = 1e-4;

dgpp = @(r) (sptm.g_uv_comp([0,r+dr,pi/2,0],4,4)-sptm.g_uv_comp([0,r-dr,pi/2,0],4,4))/(2*dr);
dgtt = @(r) (sptm.g_uv_comp([0,r+dr,pi/2,0],1,1)-sptm.g_uv_comp([0,r-dr,pi/2,0],1,1))/(2*dr);
dgtp = @(r) (sptm.g_uv_comp([0,r+dr,pi/2,0],1,4)-sptm.g_uv_comp([0,r-dr,pi/2,0],1,4))/(2*dr);

ddgpp = @(r) (dgpp(r+dr)-dgpp(r-dr))/(2*dr);
ddgtt = @(r) (dgtt(r+dr)-dgtt(r-dr))/(2*dr);
ddgtp = @(r) (dgtp(r+dr)-dgtp(r-dr))/(2*dr);

L = @(r) (sptm.g_uv_comp([0,r,pi/2,0],4,4)*disk_om(r,proret,sptm) + sptm.g_uv_comp([0,r,pi/2,0],1,4));
E = @(r)- (sptm.g_uv_comp([0,r,pi/2,0],1,4)*disk_om(r,proret,sptm) + sptm.g_uv_comp([0,r,pi/2,0],1,1));

find_r = @(r) ddgtt(r)*L(r)^2 + 2*ddgtp(r)*L(r)*E(r) + ddgpp(r) * E(r)^2 - ...
    2*(-sptm.g_uv_comp([0,r,pi/2,0],4,4)*disk_om(r,proret,sptm)^2 - 2*sptm.g_uv_comp([0,r,pi/2,0],1,4)*disk_om(r,proret,sptm) - sptm.g_uv_comp([0,r,pi/2,0],1,1));

condi = true;
while condi == true
    rr = linspace(0,2*r_guess,500);
    plfind = zeros(length(rr),1);
    for i = 1:length(rr)
        plfind(i) = find_r(rr(i));
    end

    risco_test = fzero(find_r,r_guess);
    
    plot(rr,plfind)
    hold on
    scatter(risco_test,0,'r','filled')
    plot([0 2*r_guess],[0 0],'k')
    ylim([-0.1,0.1])
    hold off

    fprintf('the r ISCO value is %.3f\n',risco_test)
    fprintf("Check the plot, and make sure that the ISCO is located at the outer root point.\n" + ...
        "If its not, then change r_guess close to the outer root point.\n")
    change = input('change r_guess? (type y to change, enter to continue) : ','s');
    if change == 'y'
        r_guess = input('input new r_guess :')
    else
        condi = false;
    end
end
risco = risco_test;