#include <cstdlib>
#include <iostream>
#include <ctime>
#define NQUOTES 32
void quote(void)
{
	time_t timer;
	srand(time(&timer));
	int nq = (int)(rand() % NQUOTES);
	std::cout << std::endl;
	switch (nq)
	{
		case 0:
		{
			std::cout << "\"Imagination is more important than knowledge.\" (A. Einstein)" << std::endl;
			break;
		}
		case 1:
		{
			std::cout << "\"Everything should be made as simple as possible, but not simpler.\" (A. Einstein)" << std::endl;
			break;
		}
		case 2:
		{
			std::cout << "\"The important thing is not to stop questioning. Curiosity has its own reason for existing.\" (A. Einstein)" << std::endl;
			break;
		}
		case 3:
		{
			std::cout << "\"I am convinced that He (God) does not play dice.\" (A. Einstein)" << std::endl;
			break;
		}
		case 4:
		{
			std::cout << "\"Peace cannot be kept by force. It can only be achieved by understanding.\" (A. Einstein)" << std::endl;
			break;
		}
		case 5:
		{
			std::cout << "\"In order to form an immaculate member of a flock of sheep one must, above all, be a sheep.\" (A. Einstein)" << std::endl;
			break;
		}
		case 6:
		{
			std::cout << "\"Anyone who has never made a mistake has never tried anything new.\" (A. Einstein)" << std::endl;
			break;
		}
		case 7:
		{
			std::cout << "\"Great spirits have often encountered violent opposition from weak minds.\" (A. Einstein)" << std::endl;
			break;
		}
		case 8:
		{
			std::cout << "\"Gravitation is not responsible for people falling in love.\" (A. Einstein)" << std::endl;
			break;
		}
		case 9:
		{
			std::cout << "\"God does not care about our mathematical difficulties. He integrates empirically.\" (A. Einstein)" << std::endl;
			break;
		}
		case 10:
		{
			std::cout << "\"Do not worry about your difficulties in Mathematics. I can assure you mine are still greater.\" (A. Einstein)" << std::endl;
			break;
		}
		case 11:
		{
			std::cout << "\"Put your hand on a hot stove for a minute, and it seems like an hour. Sit with a pretty girl for an hour, and it seems like a minute. THAT'S relativity.\" (A. Einstein)" << std::endl;
			break;
		}
		case 12:
		{
			std::cout << "\"Any body suspended in space will remain in space until made aware of its situation. Daffy Duck steps off a cliff, expecting further pastureland. He loiters in midair, soliloquizing flippantly, until he chances to look down. At this point, the familiar principle of 32 feet per second per second takes over.\" (Physics of Cartoons - Law I)" << std::endl;
			break;
		}
		case 13:
		{
			std::cout << "\"Any body in motion will tend to remain in motion until solid matter intervenes suddenly. Whether shot from a cannon or in hot pursuit on foot, cartoon characters are so absolute in their momentum that only a telephone pole or an outsize boulder retards their forward motion absolutely. Sir Isaac Newton called this sudden termination of motion the stooge's surcease.\" (Physics of Cartoons - Law II)" << std::endl;
			break;
		}
		case 14:
		{
			std::cout << "\"Any body passing through solid matter will leave a perforation conforming to its perimeter. Also called the silhouette of passage, this phenomenon is the speciality of victims of directed-pressure explosions and of reckless cowards who are so eager to escape that they exit directly through the wall of a house, leaving a cookie-cutout-perfect hole. The threat of skunks or matrimony often catalyzes this reaction.\" (Physics of Cartoons - Law III)" << std::endl;
			break;
		}
		case 15:
		{
			std::cout << "\"The time required for an object to fall twenty stories is greater than or equal to the time it takes for whoever knocked it off the ledge to spiral down twenty flights to attempt to capture it unbroken. Such an object is inevitably priceless, the attempt to capture it inevitably unsuccessful.\" (Physics of Cartoons - Law IV)" << std::endl;
			break;
		}
		case 16:
		{
			std::cout << "\"All principles of gravity are negated by fear. Psychic forces are sufficient in most bodies for a shock to propel them directly away from the earth's surface. A spooky noise or an adversary's signature sound will induce motion upward, usually to the cradle of a chandelier, a treetop, or the crest of a flagpole. The feet of a character who is running or the wheels of a speeding auto need never touch the ground, especially when in flight.\" (Physics of Cartoons - Law V)" << std::endl;
			break;
		}
		case 17:
		{
			std::cout << "\"As speed increases, objects can be in several places at once. This is particularly true of tooth-and-claw fights, in which a character's head may be glimpsed emerging from the cloud of altercation at several places simultaneously. This effect is common as well among bodies that are spinning or being throttled. A `wacky' character has the option of self- replication only at manic high speeds and may ricochet off walls to achieve the velocity required.\" (Physics of Cartoons - Law VI)" << std::endl;
			break;
		}
		case 18:
		{
			std::cout << "\"Certain bodies can pass through solid walls painted to resemble tunnel entrances; others cannot. This trompe l'oeil inconsistency has baffled generations, but at least it is known that whoever paints an entrance on a wall's surface to trick an opponent will be unable to pursue him into this theoretical space. The painter is flattened against the wall when he attempts to follow into the painting. This is ultimately a problem of art, not of science.\" (Physics of Cartoons - Law VII)" << std::endl;
			break;
		}
		case 19:
		{
			std::cout << "\"Any violent rearrangement of feline matter is impermanent. Cartoon cats possess even more deaths than the traditional nine lives might comfortably afford. They can be decimated, spliced, splayed, accordion-pleated, spindled, or disassembled, but they cannot be destroyed. After a few moments of blinking self pity, they reinflate, elongate, snap back, or solidify. Corollary: A cat will assume the shape of its container.\" (Physics of Cartoons - Law VIII)" << std::endl;
			break;
		}
		case 20:
		{
			std::cout << "\"Everything falls faster than an anvil.\"  (Physics of Cartoons - Law IX)" << std::endl;
			break;
		}
		case 21:
		{
			std::cout << "\"For every vengeance there is an equal and opposite revengeance. This is the one law of animated cartoon motion that also applies to the physical world at large. For that reason, we need the relief of watching it happen to a duck instead.\" (Physics of Cartoons - Law X)" << std::endl;
			break;
		}
		case 22:
		{
			std::cout << "\"A sharp object will always propel a character upward. When poked (usually in the buttocks) with a sharp object (usually a pin), a character will defy gravity by shooting straight up, with great velocity.\" (Physics of Cartoons - Law Amendment A)" << std::endl;
			break;
		}
		case 23:
		{
			std::cout << "\"The laws of object permanence are nullified for \"cool\" characters. Characters who are intended to be \"cool\" can make previously nonexistent objects appear from behind their backs at will. For instance, the Road Runner can materialize signs to express himself without speaking.\" (Physics of Cartoons - Law Amendment B)" << std::endl;
			break;
		}
		case 24:
		{
			std::cout << "\"Explosive weapons cannot cause fatal injuries. They merely turn characters temporarily black and smoky.\" (Physics of Cartoons - Law Amendment C)" << std::endl;
			break;
		}
		case 25:
		{
			std::cout << "\"Gravity is transmitted by slow-moving waves of large wavelengths. Their operation can be wittnessed by observing the behavior of a canine suspended over a large vertical drop. Its feet will begin to fall first, causing its legs to stretch. As the wave reaches its torso, that part will begin to fall, causing the neck to strech. As the head begins to fall, tension is released and the canine will resume its regular proportions until such time as it strikes the ground.\" (Physics of Cartoons - Law Amendment D)" << std::endl;
			break;
		}
		case 26:
		{
			std::cout << "\"Dynamite is spontaneously generated in \"C-spaces\" (spaces in which cartoon laws hold). The process is analogous to steady-state theories of the universe which postulated that the tensions involved in maintianing a space would cause the creation of hydrogen from nothing. Dynamite quanta are quite large (stick sized) and unstable (lit). Such quanta are attracted to psychic forces generated by feelings of distress in \"cool\" characters, who are able to use said quanta to their advantage. One may imagine C-spaces where all matter and energy result from primal masses of dynamite exploding. A big bang indeed.\" (Physics of Cartoons - Law Amendment E)" << std::endl;
			break;
		}
		case 27:
		{
			std::cout << "\"F = ma\" (I. Newton)" << std::endl;
			break;
		}
		case 28:
		{
			std::cout << "\"Eppur si muove (Yet it moves).\" (G. Galilei)" << std::endl;
			break;
		}
		case 29:
		{
			std::cout << "\"And all this science I can't understand it's just my job five days a week\" (\"Rocket man\" - E. John)" << std::endl;
			break;
		}
		case 30:
		{
			std::cout << "\"Why did the chicken cross the road?\" \"How many roads must one chicken cross?\" (Bob Dylan)" << std::endl;
			break;
		}
		case 31:
		{
			std::cout << "\"Grande Giove!\" (Doc E. L. Brown - \"Back to the future\")" << std::endl;
			break;
		}
	}
	std::cout << std::endl;
	return;
}
