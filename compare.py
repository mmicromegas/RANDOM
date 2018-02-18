import NETWORKS
import matplotlib.pyplot as plt

pfl = 'prompi_t2e9d1e7.dat'

aprox13_fl1 = '/home/miro/cococubed/aprox13/aprox13_t2e9d1e7_cf88_options0010111_0000_z1.dat'
aprox13_fl2 = '/home/miro/cococubed/aprox13/aprox13_t2e9d1e7_cf88_options0010111_0002_z1.dat'

#aprox13_fl1 = '/home/miro/cococubed/aprox13/aprox13_t2e9d1e7_jina_options0100111_0000_z1.dat'
#aprox13_fl2 = '/home/miro/cococubed/aprox13/aprox13_t2e9d1e7_jina_options0100111_0002_z1.dat'

torch46_fl1 = '/home/miro/cococubed/torch/torch46iso_t2e9d1e7_cf88_options01100111_0000_z1.dat'
torch46_fl2 = '/home/miro/cococubed/torch/torch46iso_t2e9d1e7_cf88_options01100111_0002_z1.dat'
torch46_fl3 = '/home/miro/cococubed/torch/torch46iso_t2e9d1e7_cf88_options01100111_0003_z1.dat'
torch46_fl4 = '/home/miro/cococubed/torch/torch46iso_t2e9d1e7_cf88_options01100111_0004_z1.dat'
torch46_fl5 = '/home/miro/cococubed/torch/torch46iso_t2e9d1e7_cf88_options01100111_0005_z1.dat'

#oburn25_fl1 = '/home/miro/cococubed/torch/test40000_z1.dat'
#oburn25_fl2 = '/home/miro/cococubed/torch/test40002_z1.dat'
#oburn25_fl3 = '/home/miro/cococubed/torch/test40003_z1.dat'
#oburn25_fl4 = '/home/miro/cococubed/torch/test40004_z1.dat'

prompi  = NETWORKS.PROMPI_BURN(pfl)

aprox13 = NETWORKS.COCOCUBED(aprox13_fl1,aprox13_fl2)

torch46 = NETWORKS.COCOCUBED_TORCH(torch46_fl1,torch46_fl2,\
                                   torch46_fl3,torch46_fl4,torch46_fl5)

#oburn25 =  NETWORKS.COCOCUBED_OBURN25(oburn25_fl1,oburn25_fl2,\
#                                      oburn25_fl3,oburn25_fl4)

timep  = prompi.data('time')
timeaprox13 = aprox13.datae('time')
timetorch46 = torch46.datae('time')
#timeoburn25 = oburn25.datae('time')

#start = time[timep.index[1]]
#stop  = time[timep.index[-1]]
        
fig, ax1 = plt.subplots(figsize=(7,6))

start = 1.e-5
stop  = 1.e-1
yb    = 1.e17
yt    = 1.e20
        
ax1.axis([start,stop,yb,yt])

ax1.loglog(timep,prompi.data('eps')-prompi.data('snu'),label='prompi oburn25 s-snu')
ax1.loglog(timeaprox13,aprox13.datae('s-snu'),label='aprox 13 cf88 s-snu')
ax1.loglog(timetorch46,torch46.datae('s-snu'),marker='+',label='torch 46 cf88 s-snu')
#ax1.loglog(timeoburn25,oburn25.datae('s-snu'),label='oburn25 s-snu')

ax1.set_xlabel('time (s)')
ax1.set_ylabel('e (erg g$^{-1}$ s$^{-1}$)')
ax1.legend(loc=4,prop={'size':8})

plt.title('X(c12) = 0.5, X(o16) = 0.5, t2e9d1e7 noscreen')
plt.show(block=False) 
plt.savefig('ozb_enuc_t2e9d1e7_noscreen_CF88.png')

fig, ax2 = plt.subplots(figsize=(7,6))

#start = 1.e-10
#stop  = 1.e-6
yb    = 1.e-3
yt    = 1.e0

ax2.axis([start,stop,yb,yt])

ax2.loglog(timep,prompi.data('he4'),label='prompi oburn25 he4',color='y')
ax2.loglog(timep,prompi.data('c12'),label='prompi oburn25 c12',color='r')
ax2.loglog(timep,prompi.data('o16'),label='prompi oburn25 o16',color='g')
ax2.loglog(timep,prompi.data('ne20'),label='prompi oburn25 ne20',color='b')
ax2.loglog(timep,prompi.data('mg24'),label='prompi oburn25 mg24',color='k')
ax2.loglog(timep,prompi.data('si28'),label='prompi oburn25 si28',color='c')

ax2.loglog(timeaprox13,aprox13.datax('he4'),label='aprox13 cf88 he4',linestyle=':',color='y')
ax2.loglog(timeaprox13,aprox13.datax('c12'),label='aprox13 cf88 c12',linestyle=':',color='r')
ax2.loglog(timeaprox13,aprox13.datax('o16'),label='aprox13 cf88 o16',linestyle=':',color='g')
ax2.loglog(timeaprox13,aprox13.datax('ne20'),label='aprox13 cf88 ne20',linestyle=':',color='b')
ax2.loglog(timeaprox13,aprox13.datax('mg24'),label='aprox13 cf88 mg24',linestyle=':',color='k')
ax2.loglog(timeaprox13,aprox13.datax('si28'),label='aprox13 cf88 si28',linestyle=':',color='c')

ax2.loglog(timetorch46,torch46.datax5('he4'),label='torch46 cf88 he4',marker='+',color='y')
ax2.loglog(timetorch46,torch46.datax2('c12'),label='torch46 cf88 c12',marker='+',color='r')
ax2.loglog(timetorch46,torch46.datax3('o16'),label='torch46 cf88 o16',marker='+',color='g')
ax2.loglog(timetorch46,torch46.datax3('ne20'),label='torch46 cf88 ne20',marker='+',color='b')
ax2.loglog(timetorch46,torch46.datax3('mg24'),label='torch46 cf88 mg24',marker='+',color='k')
ax2.loglog(timetorch46,torch46.datax4('si28'),label='torch46 cf88 si28',marker='+',color='c')


ax2.set_xlabel('time (s)')
ax2.set_ylabel('X')
ax2.legend(loc=4,prop={'size':5})

plt.title('X(c12) = 0.5, X(o16) = 0.5, t2e9d1e7 noscreen')		
plt.show(block=False)        
plt.savefig('ozb_x_t2e9d1e7_noscreen_CF88.png')
