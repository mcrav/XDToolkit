'''
#####################################################################
#-------------------EMAIL--------------------------------------------
#####################################################################
'''
import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.application import MIMEApplication
from email.mime.text import MIMEText
import os

def sendEmail(body = '', email = '', attachments = [], subject = '', toaddr = ''):
    '''
    Send email to my email account from xdtoolkit@gmail.com
    '''
    print('Sending email')
    fromaddr = "xdtoolkit@gmail.com"
    if not toaddr:
        toaddr = "mcrav@chem.au.dk"
    msg = MIMEMultipart()
    msg['From'] = fromaddr
    msg['To'] = toaddr
    msg['Subject'] = subject

    for f in attachments:

        with open(f, "rb") as fil:
            part = MIMEApplication(
                fil.read(),
                Name=os.path.basename(f)
            )
            part['Content-Disposition'] = 'attachment; filename="%s"' % os.path.basename(f)
            msg.attach(part)

    bodyText = '{0}<br><br>Email Address: {1}'.format(body,email)
    msg.attach(MIMEText(bodyText, 'html'))
    print('starting smtp')
    server = smtplib.SMTP('smtp.gmail.com', 587)
    server.starttls()
    server.login("xdtoolkit@gmail.com", '****')
    print('logged in')
    text = msg.as_string()
    server.sendmail("xdtoolkit@gmail.com", toaddr, text)
    print('sent')

    server.quit()
