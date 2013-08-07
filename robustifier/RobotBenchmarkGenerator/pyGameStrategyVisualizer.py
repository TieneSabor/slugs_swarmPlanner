#!/usr/bin/python
#
# Performs \infty-resilient GR(1) synthesis
# Takes the input file in Slugs format

import math
import os
import sys, code
import resource
import subprocess
import signal
import tempfile
import copy
import itertools
import Image
import os, pygame, pygame.locals

# ==================================
# Settings
# ==================================
MAGNIFY = 64

# ==================================
# Entry point
# ==================================
if len(sys.argv)<2:
    print >>sys.stderr, "Error: Need PNG file as parameter"
    sys.exit(1)
specFile = sys.argv[1]

# ==================================
# Read input image
# ==================================
import os,sys
pngfile = Image.open(specFile)
pngFileBasis = specFile[0:specFile.rfind(".png")]
# print "Size of Workspace:",pngfile.size
xsize = pngfile.size[0]
ysize = pngfile.size[1]
imageData = pngfile.getdata()
palette = pngfile.getpalette()
if (xsize>1023):
    print >>sys.stderr,"Error: Scenario is too large - not supported."
    sys.exit(1)
if (ysize>1023):
    print >>sys.stderr,"Error: Scenario is too large - not supported."
    sys.exit(1)

# Adapt size of display
pygame.init()
displayInfo = pygame.display.Info()
MAGNIFY = min(MAGNIFY,displayInfo.current_w*3/4/xsize)
MAGNIFY = min(MAGNIFY,displayInfo.current_h*3/4/ysize)


# ====================================================================
# Read input image color encoding configuration file
# ====================================================================
colorCodingInformationFile = open(pngFileBasis+".colorcoding","r").readlines()
colorCoding = [(int(a.strip().split(" ")[0]),a.strip().split(" ")[1]) for a in colorCodingInformationFile if a.strip()!=""]
colorCodingMap = {a:b for (a,b) in colorCoding}
print colorCoding

# ====================================================================
# Assign keys to doors and deliveries
# ====================================================================
keys = []
for (a,b) in colorCoding:
    if b=="Door":
        keys.append((a,b))
    elif b=="Delivery":
        keys.append((a,b))

# ====================================================================
# Assign robot keys
# ====================================================================
movingObstacles = []
for (a,b) in colorCoding:
    if b.startswith("MovingObstacle:"):
        movingObstacles.append((a,b.split(":")))

# ==================================
# Prepare Slugs Call
# ==================================
slugsinfile = specFile[0:specFile.rfind(".png")]+".slugsin"
print slugsinfile
slugsLink = sys.argv[0][0:sys.argv[0].rfind("pyGameStrategyVisualizer.py")]+"../../src/slugs"
print slugsLink

# ==================================
# Main loop
# ==================================
def actionLoop():
    screen = pygame.display.set_mode(((xsize+2)*MAGNIFY,(ysize+2)*MAGNIFY))
    pygame.display.set_caption('Strategy visualizer')
    clock = pygame.time.Clock()

    screenBuffer = pygame.Surface(screen.get_size())
    screenBuffer = screenBuffer.convert()
    screenBuffer.fill((64, 64, 64)) # Dark Gray

    # Open Slugs
    slugsProcess = subprocess.Popen(slugsLink+" --interactiveStrategy "+slugsinfile, shell=True, bufsize=1048000, stdin=subprocess.PIPE, stdout=subprocess.PIPE)

    # Get input APs
    slugsProcess.stdin.write("XPRINTINPUTS\n")
    slugsProcess.stdin.flush()
    slugsProcess.stdout.readline() # Skip the prompt
    lastLine = " "
    inputAPs = []
    while (lastLine!=""):
        lastLine = slugsProcess.stdout.readline().strip()
        if lastLine!="":
            inputAPs.append(lastLine)

    # Get output APs
    slugsProcess.stdin.write("XPRINTOUTPUTS\n")
    slugsProcess.stdin.flush()
    slugsProcess.stdout.readline() # Skip the prompt
    lastLine = " "
    outputAPs = []
    while (lastLine!=""):
        lastLine = slugsProcess.stdout.readline().strip()
        if lastLine!="":
            outputAPs.append(lastLine)
    print inputAPs
    print outputAPs

    # Get initial state
    slugsProcess.stdin.write("XGETINIT\n")
    slugsProcess.stdin.flush()
    slugsProcess.stdout.readline() # Skip the prompt
    currentState = slugsProcess.stdout.readline().strip()
    print "currentState: "+currentState

    # Pre-store positions
    doorAndDeliveryInputBitPositions = {}
    for (a,b) in colorCoding:
        if b=="Door" or b=="Delivery":
            for pos,name in enumerate(inputAPs):
                if name=="door"+str(a) or name=="deliveryrequest"+str(a)  :
                    doorAndDeliveryInputBitPositions[a] = pos

    while 1:

        for event in pygame.event.get():
            if event.type == pygame.locals.QUIT or (event.type == pygame.locals.KEYDOWN and event.key == pygame.locals.K_ESCAPE):
                slugsProcess.stdin.write("QUIT\n")
                slugsProcess.stdin.flush()
                return

        # Obtain robot information for drawing
        robotX = 0
        robotY = 0
        robotDeltaX = 0
        robotDeltaY = 0
        for i,ap in enumerate(inputAPs):
            if ap in ["x0","x1","x2","x3","x4","x5","x6","x7","x8","x9"]:
                if currentState[i]=="1":
                    robotX += (1 << int(ap[1:]))
                elif currentState[i]=="0":
                    pass
                else:
                    raise 123
            elif ap in ["y0","y1","y2","y3","y4","y5","y6","y7","y8","y9"]:
                if currentState[i]=="1":
                    robotY += (1 << int(ap[1:]))

        for i,ap in enumerate(outputAPs):
            if ap=="left":
                if currentState[i+len(inputAPs)]!="0":
                    robotDeltaX -= 1
            elif ap=="right":
                if currentState[i+len(inputAPs)]!="0":
                    robotDeltaX += 1
            elif ap=="up":
                if currentState[i+len(inputAPs)]!="0":
                    robotDeltaY -= 1
            elif ap=="down":
                if currentState[i+len(inputAPs)]!="0":
                    robotDeltaY += 1

        # Draw Field
        for x in xrange(0,xsize):
            for y in xrange(0,ysize):
                paletteColor = imageData[y*xsize+x]
                
                # Translate color to what is to be written
                if paletteColor in colorCodingMap:
                    colorCodeDescription = colorCodingMap[paletteColor]
                    if colorCodeDescription.startswith("MovingObstacle:"):
                        color = (255,255,255)
                    elif colorCodeDescription=="Door" or colorCodeDescription=="Delivery":
                        if currentState[doorAndDeliveryInputBitPositions[paletteColor]]=="0":
                            color = (128+palette[paletteColor*3]/2,128+palette[paletteColor*3+1]/2,128+palette[paletteColor*3+2]/2)
                        else:
                            color = palette[paletteColor*3:paletteColor*3+3]
                    else:
                        color = palette[paletteColor*3:paletteColor*3+3]
                else:
                    color = palette[paletteColor*3:paletteColor*3+3]

                pygame.draw.rect(screenBuffer,color,((x+1)*MAGNIFY,(y+1)*MAGNIFY,MAGNIFY,MAGNIFY),0)

        # Draw "Good" Robot
        pygame.draw.circle(screenBuffer, (192,32,32), ((robotX+1)*MAGNIFY+MAGNIFY/2,(robotY+1)*MAGNIFY+MAGNIFY/2) , MAGNIFY/3-2, 0)
        pygame.draw.circle(screenBuffer, (255,255,255), ((robotX+1)*MAGNIFY+MAGNIFY/2,(robotY+1)*MAGNIFY+MAGNIFY/2) , MAGNIFY/3-1, 1)
        pygame.draw.circle(screenBuffer, (0,0,0), ((robotX+1)*MAGNIFY+MAGNIFY/2,(robotY+1)*MAGNIFY+MAGNIFY/2) , MAGNIFY/3, 1)

        for x in xrange(0,xsize):
            for y in xrange(0,ysize):
                pygame.draw.rect(screenBuffer,(0,0,0),((x+1)*MAGNIFY,(y+1)*MAGNIFY,MAGNIFY,MAGNIFY),1)
        pygame.draw.rect(screenBuffer,(0,0,0),(MAGNIFY-1,MAGNIFY-1,MAGNIFY*xsize+2,MAGNIFY*ysize+2),1)

        # Flip!
        screen.blit(screenBuffer, (0, 0))
        pygame.display.flip()


        # Update Doors and requests
        nextInput = currentState[0:len(inputAPs)]
        for keyNum,(a,b) in enumerate(keys):
            if pygame.key.get_pressed()[pygame.locals.K_1+keyNum]:
                nextInput = nextInput[0:doorAndDeliveryInputBitPositions[a]]+"1"+nextInput[doorAndDeliveryInputBitPositions[a]+1:]
            else:
                nextInput = nextInput[0:doorAndDeliveryInputBitPositions[a]]+"0"+nextInput[doorAndDeliveryInputBitPositions[a]+1:]

        # Update robot X and Y positions
        robotX = robotX + robotDeltaX
        robotY = robotY + robotDeltaY
        print "Pos and acc: ",robotX, robotY, robotDeltaX, robotDeltaY
        for i,ap in enumerate(inputAPs):
            if ap in ["x0","x1","x2","x3","x4","x5","x6","x7","x8","x9"]:
                if (robotX & (1 << int(ap[1:])))>0:
                    nextInput = nextInput[0:i]+"1"+nextInput[i+1:]
                else:
                    nextInput = nextInput[0:i]+"0"+nextInput[i+1:]
            if ap in ["y0","y1","y2","y3","y4","y5","y6","y7","y8","y9"]:
                if (robotY & (1 << int(ap[1:])))>0:
                    nextInput = nextInput[0:i]+"1"+nextInput[i+1:]
                else:
                    nextInput = nextInput[0:i]+"0"+nextInput[i+1:]

        print currentState, nextInput

        # Make the transition
        slugsProcess.stdin.write("XMAKETRANS\n"+nextInput)
        slugsProcess.stdin.flush()
        print slugsProcess.stdout.readline() # Skip the prompt
        nextLine = slugsProcess.stdout.readline().strip()
        print nextLine
        if nextLine=="ERROR":
            screenBuffer.fill((192, 64, 64)) # Red!
            # Keep the state the same
        else:
            currentState = nextLine
            screenBuffer.fill((64, 64, 64)) # Gray, as usual

        # Done
        clock.tick(10)


# ==================================
# Call main program
# ==================================
actionLoop()
