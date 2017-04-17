'''
#####################################################################
#-------------------EMAIL--------------------------------------------
#####################################################################
'''

def sendEmail(body = '', email = '', attachments = [], subject = ''):
    '''
    Send email to my email account from xdtoolkit@gmail.com
    '''
    fromaddr = "xdtoolkit@gmail.com"
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

    server = smtplib.SMTP('smtp.gmail.com', 587)
    server.starttls()
    server.login("xdtoolkit@gmail.com", '***')
    text = msg.as_string()
    server.sendmail("xdtoolkit@gmail.com", "mcrav@chem.au.dk", text)

    server.quit()
