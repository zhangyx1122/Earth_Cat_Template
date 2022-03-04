package com.company;

import java.io.*;

public class Main {
    static FileWriter fw;
    static BufferedWriter bufw;
    private String line;

    public static void fileOutput(File dir) throws IOException {

        FileReader fr = new FileReader(dir.toString());
        BufferedReader bufr = new BufferedReader(fr);

        String line = null;
        while ((line = bufr.readLine()) != null) {
            bufw.write(line);
            bufw.newLine();
            bufw.flush();
        }


        bufr.close();
    }

    public static void visitAllDirsAndFiles(File dir, int num) throws IOException {
        String pre = "";
        for (int i = 0; i < num; ++i) pre = pre + "#";
        String s = dir.toString();
        s = s.substring(s.lastIndexOf('\\') + 1);
        String suff = s.substring(s.lastIndexOf(".") + 1);
        if (suff.equals("png") || suff.equals("jpg")) return;
        boolean f = false;
        if (s.indexOf('.') > 0) {
            s = s.substring(0, s.indexOf('.'));
            f = true;
        }

//        System.out.print(pre);
//        System.out.println(s);
        bufw.write(pre);
        bufw.write(s);
        bufw.newLine();
        bufw.flush();
        if (suff.equals("cpp")) {
            bufw.newLine();
            bufw.write("\n```\n");
            bufw.newLine();
            bufw.flush();
        }
        if (!dir.isDirectory() && f) fileOutput(dir);
        if (suff.equals("cpp")) {
            bufw.newLine();
            bufw.write("\n```\n");
            bufw.newLine();
            bufw.flush();
        }

        if (dir.isDirectory()) {
            String[] children = dir.list();
            for (int i = 0; i < children.length; i++) {
                visitAllDirsAndFiles(new File(dir, children[i]), num + 1);
            }
        }
    }

    public static void main(String[] args) throws IOException {

        fw = new FileWriter("C:\\Users\\61049\\Desktop\\test.md");
        bufw = new BufferedWriter(fw);

        File dir = new File("C:\\Users\\61049\\Desktop\\testFile");

        visitAllDirsAndFiles(dir, 0);

        bufw.close();
    }
}
