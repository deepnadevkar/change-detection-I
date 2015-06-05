Attribute VB_Name = "Module2"
Option Explicit

'State Variables
Global gsInitState As String

Global Const TOUCH_BASE_ADDRESS = &H3F8  '&H300
Global Const TOUCH_DATA_REGISTER = TOUCH_BASE_ADDRESS
Global Const TOUCH_STATUS_REGISTER = TOUCH_BASE_ADDRESS + 1
Global Const TOUCH_RESET_REGISTER = TOUCH_BASE_ADDRESS + 2

Global Const TOUCH_READY_MASK = &H1
Global Const TOUCH_DATA_MASK = &H2

'Communications
Global Const TOUCH_HARDWARE_FLOW_CONTROL_ON = &H41
Global Const TOUCH_HARDWARE_FLOW_CONTROL_OFF = &H42
Global Const TOUCH_REPORT_TRANSFER_OFF = &H43
Global Const TOUCH_REPORT_TRANSFER_ON = &H44
Global Const TOUCH_GET_ONE_REPORT = &H46

'Reporting
Global Const TOUCH_SCAN_REPORTING = &H22
Global Const TOUCH_COORDINATE_REPORTING = &H23

'Touch Mode
Global Const TOUCH_ENTER_POINT_MODE = &H25
Global Const TOUCH_TRACKING_MODE = &H26
Global Const TOUCH_CONTINUOUS_MODE = &H27
Global Const TOUCH_EXIT_POINT_MODE = &H28
Global Const TOUCH_ADD_EXIT_POINT_MODIFIER = &H29

'Information Requests
Global Const TOUCH_GET_ERROR_REPORT = &H32
Global Const TOUCH_GET_CONFIGURATION_REPORT = &H33
Global Const TOUCH_GET_FIRMWARE_VERSION_REPORT = &H34
Global Const TOUCH_GET_FAILED_BEAM_REPORT = &H36
Global Const TOUCH_GET_FRAME_SIZE_REPORT = &H37
Global Const TOUCH_GET_STATE_REPORT = &H47

'System
Global Const TOUCH_ECHO_ON = &H20
Global Const TOUCH_ECHO_OFF = &H21
Global Const TOUCH_SCANNING_ON = &H2A
Global Const TOUCH_SCANNING_OFF = &H2B
Global Const TOUCH_RUN_DIAGNOSTICS = &H3A
Global Const TOUCH_SOFTWARE_RESET = &H3C
Global Const TOUCH_CLEAR_TOUCH_REPORT_BUFFER = &H3D
Global Const TOUCH_RESET = &H45
Global Const TOUCH_SWITCH_TO_SFPII = &H65

'Protocol Command Headers
Global Const TOUCH_RPTHDR_ADD_EXIT_POINT_MODIFIER = &HFD        'FD X Y FF
Global Const TOUCH_RPTHDR_CONTINUOUS_MODE = &HFE                'FE X Y FF
Global Const TOUCH_RPTHDR_COORDINATE_REPORTING = &HFE           'FE X Y FF
Global Const TOUCH_RPTHDR_ENTER_POINT_MODE = &HFE               'FE X Y FF
Global Const TOUCH_RPTHDR_EXIT_POINT_MODE = &HFE                'FE X Y FF
Global Const TOUCH_RPTHDR_GET_CONFIGURATION_REPORT = &HF7       'F7 Nt T1...Tn FF
Global Const TOUCH_RPTHDR_GET_ERROR_REPORT = &HF8               'F8 Na A1...An Nb B1...Bn FF
Global Const TOUCH_RPTHDR_GET_FAILED_BEAM_REPORT = &HFA         'FA Nx X1 T1...Xn Tn Ny Y1 S1...Yn Sn FF
Global Const TOUCH_RPTHDR_GET_FIRMWARE_VERSION_REPORT = &HF6    'F6 Nv V1...Vn FF
Global Const TOUCH_RPTHDR_GET_FRAME_SIZE_REPORT = &HF4          'F4 Nx Ny FF
Global Const TOUCH_RPTHDR_GET_ONE_REPORT = &HF5                 'F5 FF or any other report
Global Const TOUCH_RPTHDR_GET_STATE_REPORT = &HF2               'F2 Nx S1...Sn FF
Global Const TOUCH_RPTHDR_RUN_DIAGNOSTICS = &HF8                'F8 Na A1...An Nb B1...Bn FF
Global Const TOUCH_RPTHDR_SCAN_REPORTING = &HFB                 'FB Nx Xa...Xn Ny Ya...Yn FF
Global Const TOUCH_RPTHDR_TRACKING_MODE = &HFE                  'FE X Y FF

'Protocol Command Trailer
Global Const TOUCH_REPORT_TRAILER = &HFF


Public Function gfnsGetErrorReport() As String

    Dim sResponse As String
    Dim iDataByte As Integer
    
    While (agInp(TOUCH_STATUS_REGISTER) And TOUCH_DATA_MASK) <> 0
        DoEvents
    Wend
    
    iDataByte = agInp(TOUCH_DATA_REGISTER)
    sResponse = Format(Hex(iDataByte), "00")
    
    If iDataByte = TOUCH_RPTHDR_GET_ERROR_REPORT Then
        
        While (agInp(TOUCH_STATUS_REGISTER) And TOUCH_DATA_MASK) <> 0
            DoEvents
        Wend
        
        iDataByte = agInp(TOUCH_DATA_REGISTER)
        sResponse = sResponse & " " & Format(Hex(iDataByte), "00")
        
        If iDataByte <> 0 Then
            MsgBox "The touchscreen did not initialize correctly."
        End If
        
        While (agInp(TOUCH_STATUS_REGISTER) And TOUCH_DATA_MASK) <> 0
            DoEvents
        Wend
        
        iDataByte = agInp(TOUCH_DATA_REGISTER)
        sResponse = sResponse & " " & Format(Hex(iDataByte), "00")
    
        gfnsGetErrorReport = sResponse
    
    Else
        gfnsGetErrorReport = "xxx"
        MsgBox "The incorrect message was read from the touchscreen."
    End If
    
End Function

Public Sub gsubDisableTouchscreen(tmrReadTouch As Control)

    agOutp TOUCH_DATA_REGISTER, TOUCH_SCANNING_OFF

    tmrReadTouch.Enabled = False

End Sub

Public Sub gsubEnableTouchscreen(tmrReadTouch As Control)

    agOutp TOUCH_DATA_REGISTER, TOUCH_SCANNING_ON

    tmrReadTouch.Interval = 100
    tmrReadTouch.Enabled = True

End Sub

Public Sub gsubGetXY(iRawX As Integer, iRawY As Integer)

    Dim iDataByte As Long
    
'    While (agInp(TOUCH_STATUS_REGISTER) And TOUCH_DATA_MASK) <> 0
 '       DoEvents
  '  Wend

  '  iDataByte = agInp(TOUCH_DATA_REGISTER)

  '  If iDataByte = TOUCH_RPTHDR_ENTER_POINT_MODE Then
        
  '      While (agInp(TOUCH_STATUS_REGISTER) And TOUCH_DATA_MASK) <> 0
  '          DoEvents
   '     Wend
        
    '    iRawX = agInp(TOUCH_DATA_REGISTER)
        
     '   While (agInp(TOUCH_STATUS_REGISTER) And TOUCH_DATA_MASK) <> 0
        '    DoEvents
      '  Wend
        
       ' iRawY = agInp(TOUCH_DATA_REGISTER)
        
     '   While (agInp(TOUCH_STATUS_REGISTER) And TOUCH_DATA_MASK) <> 0
     '       DoEvents
     '   Wend
        
      '  iDataByte = agInp(TOUCH_DATA_REGISTER)
    
       ' If iDataByte <> TOUCH_REPORT_TRAILER Then
       '     MsgBox "Protocol framing error."
       ' End If
        
    'End If

    'iRawX = iRawX * (SCREEN_WIDTH / 92)
    'iRawY = iRawY * (SCREEN_HEIGHT / 70)
    iRawX = 11
    iRawY = 21


End Sub
Public Sub gsubInitTouchscreen(tmrInit As Control, sInitState As String, tmrReadTouch As Control)

    tmrInit.Enabled = False
    
    If sInitState = "Init" Then
        agOutp TOUCH_RESET_REGISTER, 0
        tmrInit.Interval = 1000
        sInitState = "SoftwareReset"
        tmrInit.Enabled = True
        
    ElseIf sInitState = "SoftwareReset" Then
        agOutp TOUCH_DATA_REGISTER, TOUCH_SOFTWARE_RESET
        tmrInit.Interval = 1000
        sInitState = "ReportTransferOn"
        tmrInit.Enabled = True
        
    ElseIf sInitState = "ReportTransferOn" Then
        agOutp TOUCH_DATA_REGISTER, TOUCH_REPORT_TRANSFER_ON
        tmrInit.Interval = 1000
        sInitState = "GetErrorReport"
        tmrInit.Enabled = True
        
    ElseIf sInitState = "GetErrorReport" Then
        agOutp TOUCH_DATA_REGISTER, TOUCH_GET_ERROR_REPORT
        tmrInit.Interval = 1000
        sInitState = "CoordinateReporting"
        tmrInit.Enabled = True
        
    ElseIf sInitState = "CoordinateReporting" Then
        agOutp TOUCH_DATA_REGISTER, TOUCH_COORDINATE_REPORTING
        tmrInit.Interval = 1000
        sInitState = "EnterPointMode"
        tmrInit.Enabled = True
        
    ElseIf sInitState = "EnterPointMode" Then
        agOutp TOUCH_DATA_REGISTER, TOUCH_ENTER_POINT_MODE
        tmrInit.Interval = 1000
        sInitState = "Enable"
        tmrInit.Enabled = True
        
    ElseIf sInitState = "Enable" Then
        tmrInit.Interval = 1000
        sInitState = "Init"
        tmrInit.Enabled = False
        gsubEnableTouchscreen tmrReadTouch
        
    End If
    
End Sub








