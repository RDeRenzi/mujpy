print('Make sure you have either installed mujpy,\nor its requirements (e.g. in venv) and you run from the mujpy/example/ folder') 
from os import system, getcwd
path = getcwd() 

def isok(ok,fit_type):
    if ok==0: 
        print("{}: completed (check that output is meaningful!)".format(fit_type))
    else:
        print("*************************** {}: broken!".format(fit_type))
    return ok, fit_type

string_ok = '\nSuccessfully tested fit types (check all output!): '
string_ko = 'Broken fit types: '
print('* Running mgml.822.3_4.1_fit.py: A1 single run single group fit')
command = 'python3 '+path+"/mgml.822.3-4.1_fit.py"
ok, string1 = isok(system(command),'A1')
string_ok += ' '+string1+',' if not ok else ''
string_ko += ' '+string1+',' if ok else ''

print('* Running almlgml.822.3_4.1_fit.py: A1_calib single run single group fit')
command = 'python3 '+path+"/almgml.822.3-4.1_fit.py"
ok, string1 = isok(system(command),'A1_calib')
string_ok += ' '+string1+',' if not ok else ''
string_ko += ' '+string1+',' if ok else ''

print('* Running mgml.822.3_4.2-1.1_fit.py: A20 single run single multi group sequential fit')
command = 'python3 '+path+"/mgml.822.3-4.2-1.1_fit.py"
ok, string1 = isok(system(command),'A20')
string_ok += ' '+string1+',' if not ok else ''
string_ko += ' '+string1+',' if ok else ''

print('* Running almgml.822.3_4.2-1.1_fit.py: A20_calib single run single multi group sequential fit')
command = 'python3 '+path+"/almgml.822.3-4.2-1.1_fit.py"
ok, string1 = isok(system(command),'A20_calib')
string_ok += ' '+string1+',' if not ok else ''
string_ko += ' '+string1+',' if ok else ''

print('* Running mgml.822.3_4.2-1.1_fit.py: A21 single run single multi group global fit')
command = 'python3 '+path+"/mgml.822.3-4+2-1.1_fit.py"
ok, string1 = isok(system(command),'A21')
string_ok += ' '+string1+',' if not ok else ''
string_ko += ' '+string1+',' if ok else ''

print('* Running almgml.822.3_4.2-1.1_fit.py: A21 single run single multi group global fit')
command = 'python3 '+path+"/almgml.822.3-4+2-1.1_fit.py"
ok, string1 = isok(system(command),'A21_calib')
string_ok += ' '+string1+',' if not ok else ''
string_ko += ' '+string1+',' if ok else ''

print('* Running  mgml.822.834.3-4.1_fit.py: B1 multirun sequential single group fit')
command = 'python3 ' + path+'/mgml.822.834.3-4.1_fit.py'
ok, string1 = isok(system(command),'B1')
string_ok += ' '+string1 if not ok else ''
string_ko += ' '+string1 if ok else ''


if string_ko == 'Broken fit types: ':
    string_ko += 'None'
print('\n'+40*'+.'+string_ok+'\n'+string_ko)

