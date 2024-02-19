function y = Vgshcyttomito(gsh_cyt);


atp=3000;

y=388*gsh_cyt./(gsh_cyt+(90-62.5*atp./(atp+500))) + 1460*(1+atp./(atp+500))*gsh_cyt./(gsh_cyt+9000);


