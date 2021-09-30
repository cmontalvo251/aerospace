int Numplanets = 4;
int waitcounter = 0;
float x,y;
float[] R = {3,5,8,12,14};
int counter = -1;
color r[] = {0,0,255,255,0};
color g[] = {0,255,0,255,255};
color b[] = {255,0,0,0,255};
String lines[];
String strings[] = {"A","B","C","D","E","F"};

import ddf.minim.* ;

Minim minim;

AudioPlayer au_player1, au_player2 ;


void setup()
{
  size(800,800);
  stroke(0,0,0,0);
  background(0);
  lines = loadStrings("position.txt");
  minim = new Minim(this) ;
  au_player1 = minim.loadFile("C.wav") ;
  //au_player2 = minim.loadFile("title_track.wav") ;
  //au_player2.play();
}

void draw()
{
  background(0);
  frameRate(60);
  counter++;
  ReadNext(counter);
  if (counter >= lines.length-1)
  {
     counter = lines.length-2;
  }
  //ellipse(x,y,R,R);
}

void ReadNext(int ii)
{
  waitcounter++;
  float[] xvec = {0,0,0,0,0};
  float[] yvec = {0,0,0,0,0};
  float[] delr = {0,0,0,0,0};
  String[] pieces = split(lines[ii], ',');
  //println(pieces);
  for (int jj = 0;jj< Numplanets;jj++)
  {
   x = float(pieces[jj*2]);
   y = float(pieces[jj*2+1]);
   //print("Planet = " + jj + "\n"); 
   //print("x = " + x + " , y = " + y + "\n"); 
   fill(r[jj],g[jj],b[jj]);
   ellipse(x+width/2,y+height/2,R[jj],R[jj]);
   xvec[jj] = x;
   yvec[jj] = y;
  }
  for (int ll = 0;ll<Numplanets;ll++)
  {
    for (int jj = ll+1;jj<Numplanets;jj++)
    {
      delr[ll] = sqrt(pow(xvec[ll]-xvec[jj],2) + pow(yvec[ll]-yvec[jj],2));
     if ((delr[ll] < 10) && (waitcounter>10))
     {
       int idx = 0;
       for (int kk = 1;kk<=ll;kk++)
       {
         idx = idx + (Numplanets-kk);
       }
       if (ll != 0)
       {
         idx--;
       }
       int index = ll+jj-2*(ll) + idx;
        print(ll + " " + jj + " " + idx + " " + index + "\n");
        au_player1 = minim.loadFile(strings[index] + ".wav") ;
        print(strings[index] + "\n");
        au_player1.play();
        //minim.stop();
        waitcounter = 0;
     }
    }
  }
  
}
