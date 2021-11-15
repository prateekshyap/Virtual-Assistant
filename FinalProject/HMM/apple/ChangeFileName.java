import java.io.File;
public class ChangeFileName
{
	public static void main(String[] args)
	{
		String fileName = "1";
		String copy = fileName;
		for (int j = 0; j < 10; ++j)
		{
			fileName = "1";
			fileName += Integer.toString(j);
			fileName += ".txt";
			String newName = Integer.toString(j);
			newName += ".txt";
			File file = new File(newName);
			File rename = new File(fileName);
			file.renameTo(rename);
		}
	}
}